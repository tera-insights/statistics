<?
// Copyright 2017 Tera Insights, LLC. All Rights Reserved.

function get_inputs_for_bucket_score_GLA(array $t_args, array $inputs) {
    $glaInputs = [];
    $groupingInputNames = $t_args['group'];
    grokit_assert( is_array($gbyAttNames), 'Invalid value given for groups, expected an expression name or list of expression names');
    foreach( $inputs as $name => $type ) {
      $is_argument_for_bucket_score_gla = !in_array($name, $groupingInputNames);
      if($is_argument_for_bucket_score_gla) {
          $glaInputs[$name] = $type;
      }
    }
    return $glaInputs;
}

function get_grouping_attributes(array $t_args, array $inputs) {
  $groupingInputNames = $t_args['group'];
  $grouping_attributes = [];
  foreach( $inputs as $name => $type ) {
    $is_grouping_attribute = in_array($name, $groupingInputNames);
    if($is_grouping_attribute) {
      $grouping_attributes[$name] = $type;
    }
  }
  return $grouping_attributes;
}

function get_inner_gla(array $t_args, array $inputs, array $states) {
  $innerGLA = $t_args['aggregate'];
  grokit_assert( is_gla($innerGLA), 'Non-GLA specified as aggregate for GroupBy');
  $glaInputs = get_inputs_for_bucket_score_GLA($t_args, $inputs);
  $glaOutputs = ["hashed_key" => lookupType("BASE::INT")];
  return $innerGLA->apply($glaInputs, $glaOutputs, $states);
}

function Memory_Conscious_Hashing(array $t_args, array $inputs, array $outputs, array $states)
{
    $className = generate_name("MemoryConsciousHashing");
    $outputs = ["hashed_key" => lookupType("BASE::INT")];
    $glaInputs = get_inputs_for_bucket_score_GLA($t_args, $inputs, $outputs, $states);
    $groupingInputs = get_grouping_attributes($t_args, $inputs);
    $innerGLA = get_inner_gla($t_args);
    $minimumBucketScorePercentage = $t_args['minimumBucketScorePercentage'];
    $maxNumberOfBucketsProduced = $t_args['maxNumberOfBuckets'];
    $numberOfBuckets = $t_args['arraySize'];
    $numberOfSegments = $numberOfBuckets / 10;
    $isDebugMode = get_default($t_args, 'debug', 0) > 0;
    $tuplesPerFragment = 200000;
    $sys_headers  = ['math.h', 'armadillo', 'random', 'vector', 'stdexcept',
      'map', 'cstdlib', 'string', 'iostream', 'array'];
    $seed     = get_default($t_args, 'seed', rand());
    $user_headers = [];
    $lib_headers  = ['base\HashFct.h'];
    $libraries = $innerGLA->libraries();
    $extra        = [];
    $properties   = ['tuples'];
    $result_type  = ['fragment'];
?>

class <?=$className?>;

class <?=$className?> {
 public:
  using KeySet = std::tuple<<?=typed($inputs)?>>;
  using InnerGLA = <?=$innerGLA?>;
  using ScoreType = uint64_t;
  using ArrayType = std::array<<?=$innerGLA?>, <?=$numberOfBuckets?>>;
  using ArrayOfArraysType = std::vector<ArrayType>;
  using HashType = uint64_t;
  
  using FragmentedResultIterator = struct {
    ArrayType::const_iterator current, end;
    std::size_t index;
  };

  static constexpr std::size_t kFragmentSize = <?=$tuplesPerFragment?>;
  const float minimum_bucket_score_percentage = <?=$minimumBucketScorePercentage?>;
  const uint64_t max_number_of_buckets = <?=$maxNumberOfBuckets?>;
  static const constexpr HashType seed = <?=$seed?>;
  ScoreType total_score = 0;
  
  ArrayOfArraysType segments;

  // the iterators, only 2 elements if multi, many if fragment
  std::vector<ArrayType::const_iterator> result_iterators;

  <?=$className?> :
    segments(<?=numberOfSegments?>) {}

  HashType get_bucket_key(KeySet group) {
    return SpookyHash(Hash(group), seed);
  }

  bool bucket_will_survive(HashType key) {
    auto minimum_score = minimum_bucket_score_percentage * total_score;
    return aggregate_scores[key] > minimum_score;
  }

  int numberOfBucketsThatWillSurvive() {
    int survivors = 0;
    for (auto it = aggregate_scores.begin(); it != aggregate_scores.end(); it++) {
      if (bucket_will_survive(*it)) {
        survivors++;
      }
    }
    return survivors;
  }

  bool isTooMuchMemoryUsed() {
    return numberOfBucketsThatWillSurvive() > max_number_of_buckets;
  }

  KeySet get_group_key(<?=const_typed_ref_args($groupingInputs)?>) {
    return KeySet(<?=args($groupingInputs)?>);
  }

  HashType get_bucket_key(<?=const_typed_ref_args($groupingInputs)?>) {
    auto group = get_group_key(<?=args($groupingInputs)?>);
    return get_bucket_key(group);
  }

 public:
  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    auto index = get_bucket_key(<?=args($groupingInputs)?>);
    GLAs[index].AddItem(<?=args($glaInputs)?>);
  }

  // Merges `this->aggregate_scores` and `other->aggregate_scores`.
  void AddState(<?=$className?>& other)  {
    auto other_map = other.GetAggregateScores();
    for (auto it = GLAs.begin(); it != GLAs.end(); it++) {
      it->AddState(other_map[i]);
    }
  }

  // Construct the results
  void Finalize() {
    std::size_t num_fragment = 0;
    auto it = aggregate_scores.cbegin();
    for (std::size_t index = 0; it != aggregate_scores.cend(); ++it, index++) {
      //total_score += it->GetNextResult();
      if (index % kFragmentSize == 0) {
        result_iterators.push_back(it);
      }
    }
    result_iterators.push_back(it);
<?  if ($isDebugMode) { ?>
    std::cout << "There are " << result_iterators.size() << " iterators." << std::endl;
<?  } ?>
  }

  const ArrayType &GetAggregateScores() {
    return aggregate_scores;
  }

  int GetNumFragments() {    
    return result_iterators.size() - 1;
  }

  FragmentedResultIterator* Finalize(int fragment) const {
    auto current = result_iterators[fragment];
    auto end = result_iterators[fragment + 1];
    return new FragmentedResultIterator{current, end, 0};
  }

  bool GetNextResult(FragmentedResultIterator* it, <?=typed_ref_args($outputs)?>) {
    if (it->current == it->end) {
      return false;
    }
<?  foreach (array_keys($outputs) as $index => $name) { ?>
      <?=$name?> = std::get<<?=$index?>>(it->current->first);
<?  } ?>
    ++it->current;
    return true;
  }
};

using <?=$className?>_Iterator = <?=$className?>::FragmentedResultIterator;

<?
    return [
        'kind'           => 'GLA',
        'name'           => $className,
        'system_headers' => $sys_headers,
        'user_headers'   => $user_headers,
        'lib_headers'    => $lib_headers,
        'libraries'      => $libraries,
        'extra'          => $extra,
        'iterable'       => false,
        'input'          => $inputs,
        'output'         => $outputs,
        'result_type'    => $result_type,
    ];
}
?>