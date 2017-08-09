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
    $totalScoreMultiplier = $t_args['totalScoreMultiplier'];
    $maxNumberOfBucketsProduced = $t_args['maxNumberOfBuckets'];
    $numberOfBuckets = $t_args['arraySize'];
    $numberOfSegments = 10;
    $bucketsPerSegment = ceil($numberOfBuckets / $numberOfBuckets);
    $isDebugMode = get_default($t_args, 'debug', 0) > 0;
    $tuplesPerFragment = 200000;
    $sys_headers  = ['math.h', 'armadillo', 'random', 'vector', 'stdexcept',
      'map', 'cstdlib', 'string', 'iostream', 'array', 'list'];
    $seed     = get_default($t_args, 'seed', rand());
    $user_headers = [];
    $lib_headers  = ['base\HashFct.h', 'stats\MemoryEstimators.h'];
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
  using ArrayType = std::array<<?=$innerGLA?>, <?=$bucketsPerSegment?>>;
  using ArrayOfArraysType = std::list<ArrayType>;
  using ScoreArrayType = std::array<ScoreType, <?=$bucketsPerSegment?>>;
  using HashType = uint64_t;
  
  using FragmentedResultIterator = struct {
    ArrayType::const_iterator current, end;
    std::size_t index;
  };

  static constexpr std::size_t kFragmentSize = <?=$tuplesPerFragment?>;
  const float total_score_multiplier = <?=$totalScoreMultiplier?>;
  const uint64_t max_number_of_buckets = <?=$maxNumberOfBuckets?>;
  static const constexpr HashType seed = <?=$seed?>;
  
  ArrayOfArraysType segments;

  // the iterators, only 2 elements if multi, many if fragment
  std::vector<ArrayType::const_iterator> result_iterators;

  <?=$className?>() {
    for (size_t i = 0; i < <?=$numberOfSegments?>; i++) {
      segments.push_back(ArrayType());
    }
  }

  HashType get_bucket_key(KeySet group) {
    return SpookyHash(Hash(group), seed);
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

  std::vector<ScoreArrayType> calculate_segmented_scores() {
    std::vector<ScoreArrayType> segmented_scores(<?=$numberOfSegments?>);
    auto buckets_seen = 0;
    for (auto segment_it = segments.begin(); segment_it != segments.end(); segment_it++) {
      ScoreArrayType array_of_scores;
      for (size_t i = 0; i < segment_it->second().size() && buckets_seen < <?$numberOfBuckets?>; i++) {
        ScoreType score = segment_it->second().at(gla_it).GetNextResult();
        array_of_scores[i] = score;
        
        buckets_seen++;
      }
      segmented_scores.push_back(array_of_scores);
      segments.erase(segments.begin());
    }
    return segmented_scores;
  }

  // Construct the results
  void Finalize() {
    auto num_fragment = 0;
    auto scores = calculate_segmented_scores();
    auto total_score = get_total_score(segmented_scores, <?=$numberOfBuckets?>);
    for (auto it = scores.begin(); it != scores.end(); it++) {
      keep_if_big_enough(scores->second(), total_score, total_score_multiplier);
    }

    auto num_buckets = get_number_of_buckets(scores);
    auto buckets_seen = 0;
    for (auto it = scores.begin(); it != scores.end(); it++) {
      auto score_it = it->second().begin();
      for (auto index = 0; index < score_it->size() && buckets_seen < num_buckets_; it++) {
        if (buckets_seen % kFragmentSize == 0) {
          result_iterators.push_back(it);
        }
      }
      result_iterators.push_back(it);
    }
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