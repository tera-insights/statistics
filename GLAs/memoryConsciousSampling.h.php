<?
function Memory_Conscious_Sampling_Constant_State(array $t_args) {
    $className = $t_args['className'];
    $sys_headers  = ['Random.h',
      'boost/thread/locks.hpp',
      'boost/thread/shared_mutex.hpp',
      'atomic',
      'unordered_map'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = [];
    $minimumGroupSize = $t_args['minimumGroupSize'];
    $maximumGroupsAllowed = $t_args['maximumGroupsAllowed'];
    $initialSamplingRate = $t_args['initialSamplingRate'];
    $reductionRate = $t_args['reductionRate'];
?>

class <?=$className?>ConstantState {
 private:
  using HashType = uint64_t;
  using Map = std::unordered_map<HashType, int>;

  double global_samplingRate;
  const double reductionRate = <?=$reductionRate?>;
  const int minimum_group_size = <?=$minimumGroupSize?>;
  const int maximum_groups_allowed = <?=$maximumGroupsAllowed?>;

  // Maps hashed group key to the (effective) number of tuples we have seen for
  // this group. I use the word "effective" because this algorithm resamples if
  // it sees too many groups that will survive.
  Map frequency_map;

  boost::shared_mutex mutex;

  bool isTooMuchMemoryUsed() {
    return frequency_map.size() > maximum_groups_allowed;
  }

  bool isBernoulliSuccess(double successRate) {
    return RandDouble() < successRate;
  }

  int getNumberSampled(int populationSize, double samplingRate) {
    int sampled = 0;
    for (int trial = 0; trial < populationSize; trial++) {
      if (isBernoulliSuccess(samplingRate)) {
        sampled++;
      }
    }
    return sampled;
  }

  void resampleWithReductionRate() {
    std::vector<HashType> keys_to_erase;
    for (auto it = frequency_map.begin(); it != frequency_map.end(); it++) {
      auto key = it->first;
      auto originalCount = it->second;
      if (originalCount == 0) {
        continue;
      }
      int sampleSize = getNumberSampled(originalCount, reductionRate);
      frequency_map[key] = sampleSize;
      if (sampleSize == 0) {
        keys_to_erase.push_back(key);
      }
    }
    for (HashType key : keys_to_erase) {
      frequency_map.erase(key);
    }
  }

  void resampleMap(double samplingRate) {
    global_samplingRate = samplingRate;
    resampleWithReductionRate();
  }

  bool groupWillSurvive(HashType hash) {
    return frequency_map[hash] > minimum_group_size;
  }

  bool isKeyNew(HashType hash) const {
    return frequency_map.find(hash) == frequency_map.end();
  }

  void queue_resample_if_necessary(double rate) {
    if (isTooMuchMemoryUsed()) {
      double newSamplingRate = rate * reductionRate;
      resampleMap(newSamplingRate);
    }
  }

 public:
    friend class <?=$className?>;

    <?=$className?>ConstantState()
      : global_samplingRate(<?=$initialSamplingRate?>),
        frequency_map{} {
    }

    void AddMap(const Map &other_map) {
      boost::unique_lock<boost::shared_mutex> lock(mutex);
      double sampling_rate_during_step = global_samplingRate;
      for (auto it = other_map.begin(); it != other_map.end(); it++) {
        bool isItemSampled = isBernoulliSuccess(sampling_rate_during_step);
        if (!isItemSampled) {
          continue;
        }
        frequency_map[hash] += frequency;
      }

      queue_resample_if_necessary(sampling_rate_during_step);
    }

    bool IsGroupSurvivor(HashType hash) {
      return !isKeyNew(hash) && frequency_map[hash] > 0;
    }

    double GetSamplingRate() const {
      return global_samplingRate;
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

// This version of MCS only limits the number of produced groups.
function Memory_Conscious_Sampling(array $t_args, array $inputs, array $outputs)
{
    $className = generate_name("MemoryConsciousSampling");
    $outputs = ["hashed_key" => lookupType("BASE::INT")];

    $isDebugMode = get_default($t_args, 'debug', 0);

    $sys_headers  = ['math.h', 'armadillo', 'random', 'vector', 'stdexcept',
      'map', 'cstdlib', 'string', 'iostream', 'unordered_map'];
    $user_headers = [];
    $lib_headers  = ['base\HashFct.h'];
    $libraries    = ['armadillo', 'boost_thread'];
    $extra        = [];
    $properties   = ['tuples'];
    $result_type  = ['state'];
?>

class <?=$className?>;

<?  $state = lookupResource(
  "statistics::Memory_Conscious_Sampling_Constant_State", [
    'className' => $className,
    'minimumGroupSize' => $t_args['minimumGroupSize'],
    'maximumGroupsAllowed' => $t_args['maximumGroupsAllowed'],
    'initialSamplingRate' => $t_args['initialSamplingRate'],
    'reductionRate' => $t_args['reductionRate'],
  ]
); ?>

class <?=$className?> {
 public:
  using HashType = uint64_t;
  using Map = std::unordered_map<HashType, int>;
  using State = <?=$state?>;
  State& state;
  
  Map frequency_map;

  <?=$className?>(const State& const_state)
    :state (*(State *)(&const_state)) {
  }

  HashType chainedHash(<?=const_typed_ref_args($inputs)?>) const {
    HashType result = 0;
<?  foreach ($inputs as $index => $value) {?>
    result = SpookyHash(Hash(<?=$index?>), result);
<?}?>
    return result;
  }

  void FinalizeState() {
  }

 public:
  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    const HashType hash = chainedHash(<?=args($inputs)?>);
    frequency_map[hash]++;
  }

  void AddState(<?=$className?>& other)  {}

  bool IsGroupSurvivor(<?=const_typed_ref_args($inputs)?>) const {
    auto hash = chainedHash(<?=args($inputs)?>);
    return state.IsGroupSurvivor(hash);
  }

  void ChunkBoundary(void) {
    state.AddMap(frequency_map);
    frequency_map.clear();
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
        'result_type'    => $result_type,
        'finalize_as_state' => true,
        'generated_state' => $state,
    ];
}
?>
