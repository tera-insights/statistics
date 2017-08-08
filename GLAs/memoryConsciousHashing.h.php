<?
// Copyright 2017 Tera Insights, LLC. All Rights Reserved.

function get_inputs_for_bucket_score_GLA(array $t_args, array $inputs, array $outputs, array $states) {
    $glaInputAtts = [];
    $glaOutputAtts = [];
    $gbyAttNames = $t_args['group'];

    grokit_assert( is_array($gbyAttNames), 'Invalid value given for groups, expected an expression name or list of expression names');
    foreach( $inputs as $name => $type ) {
      $is_argument_for_bucket_score_gla = !in_array($name, $gbyAttNames);
      if($is_argument_for_bucket_score_gla) {
          $glaInputAtts[$name] = $type;
      }
    }
    foreach( $outputs as $name => $type ) {
        if( !in_array($name, $gbyAttNames) ) {
            $glaOutputAtts[$name] = $type;
        }
    }

    $innerGLA = $t_args['aggregate'];
    $innerGLA = $innerGLA->apply($glaInputAtts, $glaOutputAtts, $states);
    $libraries = $innerGLA->libraries();

    $innerRes = get_first_value( $innerGLA->result_type(), [ 'multi', 'single', 'state' ] );

    if( $innerRes == 'state' ) {
        // If the result type is state, the only output is a state object
        // containing the GLA.
        $outputName = array_keys($glaOutputAtts)[0];
        $innerOutputs = [ $outputName => lookupType('base::STATE', [ 'type' => $innerGLA ] ) ];
    }
    else {
        $innerOutputs = $innerGLA->output();

        grokit_assert(\count($innerOutputs) == \count($glaOutputAtts),
            'Expected ' . \count($glaOutputAtts) . ' outputs fromm Inner GLA, got ' .
            \count($innerOutputs));
    }

    return $glaInputAtts;
}

function get_grouping_attributes(array $t_args, array $inputs) {
  $gbyAttNames = $t_args['group'];
  $grouping_attributes = [];
  foreach( $inputs as $name => $type ) {
    $is_grouping_attribute = in_array($name, $gbyAttNames);
    if($is_grouping_attribute) {
      $grouping_attributes[$name] = $type;
    }
  }
  return $grouping_attributes;
}

function Memory_Conscious_Hashing(array $t_args, array $inputs, array $outputs, array $states)
{
    $className = generate_name("MemoryConsciousHashing");
    $outputs = ["hashed_key" => lookupType("BASE::INT")];
    $glaInputs = get_inputs_for_bucket_score_GLA($t_args, $inputs, $outputs, $states);
    $groupingInputs = get_grouping_attributes($t_args, $inputs);
    $minimumBucketScorePercentage = $t_args['minimumBucketScorePercentage'];
    $maxNumberOfBuckets = $t_args['maxNumberOfBuckets'];
    $arraySize = $t_args['arraySize'];
    $isDebugMode = get_default($t_args, 'debug', 0) > 0;
    $tuplesPerFragment = 200000;
    $sys_headers  = ['math.h', 'armadillo', 'random', 'vector', 'stdexcept',
      'map', 'cstdlib', 'string', 'iostream', 'array'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $extra        = [];
    $properties   = ['tuples'];
    $result_type  = ['fragment'];
?>

class <?=$className?>;

class <?=$className?> {
 public:
  using KeySet = std::tuple<<?=typed($inputs)?>>;
  using ScoreType = uint64_t;
  using ArrayType = std::array<uint64_t, <?=$arraySize?>>;
  using HashType = uint64_t;
  
  using FragmentedResultIterator = struct {
    ArrayType::const_iterator current, end;
    std::size_t index;
  };

  static constexpr std::size_t kFragmentSize = <?=$tuplesPerFragment?>;
  const float minimum_bucket_score_percentage = <?=$minimumBucketScorePercentage?>;
  const uint64_t max_number_of_buckets = <?=$maxNumberOfBuckets?>;

  ArrayType aggregate_scores;

  std::vector<ArrayType::const_iterator> iterators_for_fragmented_result;

  HashType get_bucket_key(KeySet group) {
    return SpookyHash(Hash(group), offset);
  }

  bool bucket_will_survive(HashType key) {
    float minimum_score = minimum_bucket_score_percentage * total_score;
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
    return numberOfBucketsThatWillSurvive() > maxNumberOfBuckets;
  }

  bool isKeyNew(KeySet key) {
    return frequency_map.find(key) == frequency_map.end();
  }

  KeySet get_group_key(<?=const_typed_ref_args($groupingInputs)?>) {
    return KeySet(<?=args($groupingInputs)?>);
  }

  ScoreType get_aggregate_score(<?=const_typed_ref_args($glaInputs)?>) {
    
  }

 public:
  // Constructor, nothing of note.
  <?=$className?>()
      : bucket_scores() {
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    const KeySet key = get_group_key(<?=args($groupingInputs)?>);
    frequency_map[key] += get_aggregate_score(<?=args($glaInputs)?>);
  }

  // Merges `this->frequency_map` and `other->frequency_map`.
  void AddState(<?=$className?>& other)  {
    auto other_map = other.GetFrequencyMap();
    for (auto it = other_map.begin(); it != other_map.end(); it++) {
      auto key = it->first;
      if (isKeyNew(key)) {
        frequency_map[key] = 0;
      }
      frequency_map[key] += it->second;
    }
  }

  // Finalize is a dummy function for this GLA and only exists because it is
  // required for a GLA of type multi. Setting return_counter is only a safety
  // check as it should already be 0.
  void Finalize() {
    return_counter = 0;
  }

  const std::map<KeySet, int> &GetFrequencyMap() {
    return frequency_map;
  }

  // The groups are traversed in order and the boundaries for each are stored in
  // the iterators vector. The number of groups per fragment is determined by
  // `$fragmentSize`. 
  int GetNumFragments() {
    std::size_t num_fragment = 0;
    auto it = frequency_map.cbegin();
    // The iterator is incremented kFragmentSize times between each boundary.
    // Note that this loop immediately pushes back an iterator pointing to the
    // minimal grouping.
    for (std::size_t index = 0; it != frequency_map.cend(); ++it, index++) {
      if (index % kFragmentSize == 0) {
        iterators.push_back(it);
      }
    }

    // This pushes back the past-the-end iterator, the upper boundary for the
    // last fragment, which can have fewer than kFragmentSize groups.
    iterators.push_back(it);
<?  if ($isDebugMode) { ?>
    std::cout << "There are " << iterators.size() << " iterators." << std::endl;
<?  } ?>
    return iterators.size() - 1;
  }

  Iterator* Finalize(int fragment) const {
    return new Iterator{iterators[fragment], iterators[fragment + 1], 0};
  }

  bool GetNextResult(Iterator* it, <?=typed_ref_args($outputs)?>) {
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

using <?=$className?>_Iterator = <?=$className?>::Iterator;

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