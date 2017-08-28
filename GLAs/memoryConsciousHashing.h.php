<?
// Copyright 2017 Tera Insights, LLC. All Rights Reserved.
function Memory_Conscious_Hashing_Constant_State(array $t_args) {
    $className = $t_args['className'];
    $sys_headers  = ['atomic', 'array', 'unordered_map', 'cstdint', 'algorithm'];
    $user_headers = [];
    $lib_headers  = ['statistics\MemoryEstimators.h'];
    $libraries    = [];
    $numberOfBuckets = $t_args['numberOfBuckets'];
    $scoreType = $t_args['scoreType'];
?>

class <?=$className?>ConstantState {
 private:
  using HashType = uint64_t;
  using ScoreType = <?=$scoreType?>;
  using AtomicScoreType = std::atomic<ScoreType>;
  using ScoreArray = std::array<ScoreType, <?=$numberOfBuckets?>>;
  using RegisterType = uint8_t;
  using RegisterArray = std::array<RegisterType, <?=$numberOfBuckets?>>;

  ScoreArray scores;
  double minimum_score;

  AtomicScoreType total_score;

  RegisterArray registers;

  double minimum_total_score_multiplier = <?=$t_args['minimumTotalScoreMultiplier']?>;

  static constexpr int bits_used_for_hashing = 35;
  static constexpr uint64_t maximum_groups_produced = <?=$t_args['maximumGroupsProduced']?>;

  uint64_t get_bucket_index(HashType hash) const {
    HashType shifted = hash >> (64 - bits_used_for_hashing);
    double ratio = (shifted * 1.0) / (1UL << bits_used_for_hashing);
    return ratio * scores.size();
  }

  uint8_t get_input_for_hyper_log_log_register(HashType hash) const {
    HashType shifted = hash << bits_used_for_hashing;
    return get_index_of_leftmost_one(shifted);
  }

 public:
    friend class <?=$className?>;

    <?=$className?>ConstantState()
      : minimum_score(0.0),
        total_score(0) {
    }

    bool IsGroupSurvivor(HashType hash) const {
      auto bucket_index = get_bucket_index(hash);
      return scores.at(bucket_index) >= minimum_score;
    }

    template <typename T>
    void AddMap(const T &map) {
      for (auto it = map.begin(); it != map.end(); it++) {
        auto bucket_index = get_bucket_index(it->first);
        long score;
        it->second.GetResult(score);
        scores[bucket_index] += score;
        registers[bucket_index] = std::max(registers[bucket_index],
          get_input_for_hyper_log_log_register(it->first));
        total_score += score;
      }
    }

    double GetMinimumScore() const {
      return minimum_score;
    }

    std::vector<RegisterType> findRegistersWithMinimumScore(ScoreType min) const {
      std::vector<RegisterType> corresponding_registers;
      for (int i = 0; i < scores.size(); i++) {
        bool isBig = scores[i] >= min;
        if (isBig) {
          corresponding_registers.push_back(registers[i]);
        }
      }
      return corresponding_registers;
    }

    uint64_t estimateCardinality(ScoreType min) const {
      auto big_registers = findRegistersWithMinimumScore(min);
      auto size = big_registers.size();
      double alpha = .7213 / (1 + 1.079 / size);
      uint64_t raw_estimate = alpha * size * size * calculate_indicator_function(big_registers);
      uint64_t small_size_upper_bound = 2.5 * size;
      uint64_t intermediate_size_upper_bound = get_upper(size) * get_correction_factor(size);
      if (raw_estimate <= small_size_upper_bound) {
        auto zeroCount = count_if(big_registers.begin(), big_registers.end(),
          [](RegisterType rt) { return rt == 0; });
        if (zeroCount > 0) {
          return size * log(size * 1.0 / zeroCount);
        } else {
          return raw_estimate;
        }
      } else if (raw_estimate <= intermediate_size_upper_bound) {
        return raw_estimate;
      } else {
        return -1 * get_upper(size) * log(1 - raw_estimate / get_upper(size));
      }
    }

    void Finalize() {
      minimum_score = total_score * minimum_total_score_multiplier;
      uint64_t current_cardinality = estimateCardinality(minimum_score);
      while (current_cardinality > maximum_groups_produced) {
        std::printf("cardinality = %ld\n", current_cardinality);
        minimum_score *= 1.1;
        current_cardinality = estimateCardinality(minimum_score);
      }
      std::printf("cardinality = %ld\n", current_cardinality);
    }
};

<?
    return [
        'kind'           => 'RESOURCE',
        'name'           => $className . 'ConstantState',
        'system_headers' => $sys_headers,
        'user_headers'   => $user_headers,
        'lib_headers'    => $lib_headers,
        'libraries'      => $libraries,
        'configurable'   => false,
    ];
}

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
    $glaInputs = get_inputs_for_bucket_score_GLA($t_args, $inputs);
    $groupingInputs = get_grouping_attributes($t_args, $inputs);
    $innerGLA = get_inner_gla($t_args, $inputs, $states);
    $sys_headers  = ['math.h', 'armadillo', 'random', 'vector', 'stdexcept',
      'map', 'cstdlib', 'string', 'iostream', 'array', 'list', 'iterator', 'limits'];
    $user_headers = [];
    $lib_headers  = ['base\HashFct.h', 'statistics\MemoryEstimators.h'];
    $libraries = $innerGLA->libraries();
    $extra        = [];
    $properties   = ['tuples'];
    $result_type  = ['fragment'];
?>

class <?=$className?>;

<?  $state = lookupResource(
  "statistics::Memory_Conscious_Hashing_Constant_State", [
    'className' => $className,
    'minimumTotalScoreMultiplier' => $t_args['minimumTotalScoreMultiplier'],
    'maximumGroupsProduced' => $t_args['maximumGroupsProduced'],
    'numberOfBuckets' => $t_args['numberOfBuckets'],
    'scoreType' => 'uint16_t',
  ]
); ?>

class <?=$className?> {
 public:
  using InnerGLA = <?=$innerGLA?>;
  using HashType = uint64_t;
  using AggregateMap = std::unordered_map<HashType, InnerGLA>;
  using State = <?=$state?>;

  static const constexpr HashType seed = 0;
  State& state;
  AggregateMap chunk_aggregate_map;

  <?=$className?>(const State& const_state)
    :state (*(State *)(&const_state)) {
  }

  bool isHashNew(HashType hash) {
    return chunk_aggregate_map.find(hash) == chunk_aggregate_map.end();
  }

 public:
  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    auto hash = chainedHash(<?=args($groupingInputs)?>);
    if (isHashNew(hash)) {
      chunk_aggregate_map[hash] = InnerGLA();
    }
    chunk_aggregate_map[hash].AddItem(<?=args($glaInputs)?>);
  }

  void AddState(<?=$className?>& other) {
    auto other_map = other.GetAggregateMap();
    for (auto it = other_map.begin(); it != other_map.end(); it++) {
      if (isHashNew(it->first)) {
        chunk_aggregate_map[it->first] = InnerGLA();
      }
      chunk_aggregate_map[it->first].AddState(it->second);
    }
  }

  const AggregateMap &GetAggregateMap() const {
    return chunk_aggregate_map;
  }

  void ChunkBoundary() {
    state.AddMap(chunk_aggregate_map);
    chunk_aggregate_map.clear();
  }

  HashType chainedHash(<?=const_typed_ref_args($groupingInputs)?>) const {
    HashType result = seed;
<?  foreach ($groupingInputs as $index => $value) {?>
    result = SpookyHash(Hash(<?=$index?>), result);
<?}?>
    return result;
  }

  bool IsGroupSurvivor(<?=const_typed_ref_args($groupingInputs)?>) const {
    auto hash = chainedHash(<?=args($groupingInputs)?>);
    return state.IsGroupSurvivor(hash);
  }

  void FinalizeState() {
    state.Finalize();
  }
};

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
        'result_type'    => ['state'],
        'finalize_as_state' => true,
        'chunk_boundary' => true,
        'generated_state' => $state,
    ];
}
?>