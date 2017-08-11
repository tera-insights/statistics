<?
// Copyright 2017 Tera Insights, LLC. All Rights Reserved.

function get_inputs_for_bucket_score_GLA(array $t_args, array $inputs) {
    $glaInputs = [];
    $groupingInputNames = $t_args['group'];
    grokit_assert( is_array($groupingInputNames), 'Invalid value given for groups, expected an expression name or list of expression names');
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
    $innerGLA = get_inner_gla($t_args, $inputs, $states);
    $minimumTotalScoreMultiplier = $t_args['minimumTotalScoreMultiplier'];
    $maxNumberOfBucketsProduced = $t_args['maxNumberOfBucketsProduced'];
    $initialNumberOfBuckets = $t_args['initialNumberOfBuckets'];
    $numberOfSegments = 10;
    $bucketsPerSegment = ceil($initialNumberOfBuckets / $numberOfSegments);
    $isDebugMode = get_default($t_args, 'debug', 0) > 0;
    $tuplesPerFragment = 200000;
    $sys_headers  = ['math.h', 'armadillo', 'random', 'vector', 'stdexcept',
      'map', 'cstdlib', 'string', 'iostream', 'array', 'list', 'iterator', 'limits'];
    $seed     = get_default($t_args, 'seed', rand());
    $user_headers = [];
    $lib_headers  = ['base\HashFct.h', 'statistics\MemoryEstimators.h'];
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
  using ScoreType = unsigned long;
  using ArrayType = std::array<<?=$innerGLA?>, <?=$bucketsPerSegment?>>;
  using ArrayOfArraysType = std::list<ArrayType>;
  using ScoreArrayType = std::array<ScoreType, <?=$bucketsPerSegment?>>;
  using HashType = uint64_t;
  using HalfHashType = uint32_t;

  static constexpr std::size_t kFragmentSize = <?=$tuplesPerFragment?>;
  const float min_total_score_multiplier = <?=$minimumTotalScoreMultiplier?>;
  const uint64_t max_number_of_buckets_produced = <?=$maxNumberOfBucketsProduced?>;
  static const constexpr HashType seed = <?=$seed?>;
  const HalfHashType max_value_for_half_of_hash = std::numeric_limits<HalfHashType>::max();
  static const uint64_t initialNumberOfBuckets = <?=$initialNumberOfBuckets?>;
  
  ArrayOfArraysType segments;

  std::vector<FragmentedResultIterator<HashType>> result_iterators;

  <?=$className?>() {
    for (size_t i = 0; i < <?=$numberOfSegments?>; i++) {
      segments.push_back(ArrayType());
    }
  }

  std::vector<ScoreArrayType> calculate_segmented_scores() {
    std::vector<ScoreArrayType> segmented_scores(<?=$numberOfSegments?>);
    auto buckets_seen = 0;
    for (auto segment_it = segments.begin(); segment_it != segments.end(); segment_it++) {
      ScoreArrayType array_of_scores;
      for (size_t i = 0; i < segment_it->size() && buckets_seen < initialNumberOfBuckets; i++) {
        long score;
        segment_it->at(i).GetResult(score);
        array_of_scores[i] = score;
        buckets_seen++;
      }
      segmented_scores.push_back(array_of_scores);
      segments.erase(segments.begin());
    }
    return segmented_scores;
  }

  int get_segment_index(HashType hash) {
    HalfHashType upper_half_of_hash = hash >> 32;
    int divisor = max_value_for_half_of_hash / (<?=$numberOfSegments?> - 1);
    return upper_half_of_hash / divisor;
  }

  int get_GLA_index(HashType hash) {
    HalfHashType lower_half_of_hash = hash << 32;
    int divisor = max_value_for_half_of_hash / (<?=$bucketsPerSegment?> - 1);
    return lower_half_of_hash / divisor;
  }

 public:
  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    auto group = KeySet(<?=args($groupingInputs)?>);
    auto hash = SpookyHash(Hash(group), seed);
    auto segment_it = segments.begin();
    std::advance(segment_it, get_segment_index(hash));
    segment_it->at(get_GLA_index(hash)).AddItem(<?=args($glaInputs)?>);
  }

  void AddState(const <?=$className?>& other)  {
    ArrayOfArraysType other_map = other.GetAggregateScores();
    for (size_t i = 0; i < other_map.size(); i++) {
      for (auto other_segment_it = other_map.begin();
        other_segment_it != other_map.end(); other_segment_it++) {
        for (size_t j = 0; j < other_segment_it->size(); j++) {
          auto beginning_it = segments.begin();
          std::advance(beginning_it, i);
          const auto my_segment_it = beginning_it;
          my_segment_it->at(j).AddState(other_segment_it->at(j));
        }
      }
    }
  }

  void Finalize() {
    auto scores = calculate_segmented_scores();
    auto total_score = get_total_score(scores, initialNumberOfBuckets);
    unsigned long minimum_score = min_total_score_multiplier * total_score;
    auto filtered = keep_if_big_enough(scores, minimum_score);
    result_iterators = build_result_iterators<FragmentedResultIterator<HashType>>(filtered, kFragmentSize);
  }

  const ArrayOfArraysType &GetAggregateScores() const {
    return segments;
  }

  int GetNumFragments() {
    return result_iterators.size() - 1;
  }

  FragmentedResultIterator<HashType>* Finalize(int fragment) {
    return &result_iterators.at(fragment);
  }

  bool GetNextResult(FragmentedResultIterator<ScoreType>* it, <?=typed_ref_args($outputs)?>) {
    if (it->current == it->end) {
      return false;
    }
    hashed_key = *(it->current);
    ++it->current;
    return true;
  }
};

using <?=$className?>_Iterator = FragmentedResultIterator<<?=$className?>::ScoreType>;

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