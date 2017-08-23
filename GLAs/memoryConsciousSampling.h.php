<?
function Memory_Conscious_Sampling_Constant_State(array $t_args) {
    $className = $t_args['className'];
    $sys_headers  = ['pthread.h', 'atomic'];
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
  using Map = std::map<HashType, int>;

  double global_samplingRate = <?=$initialSamplingRate?>;
  const double reductionRate = <?=$reductionRate?>;
  const int minimum_group_size = <?=$minimumGroupSize?>;
  const int maximum_groups_allowed = <?=$maximumGroupsAllowed?>;

  // Maps hashed group key to the (effective) number of tuples we have seen for
  // this group. I use the word "effective" because this algorithm resamples if
  // it sees too many groups that will survive.
  Map frequency_map;

  std::atomic_ulong surviving_groups;

  pthread_rwlock_t lock;

  bool isTooMuchMemoryUsed() {
    return surviving_groups > maximum_groups_allowed;
  }

  bool isBernoulliSuccess(double successRate) {
    double random = rand() * 1.0 / RAND_MAX;
    return random < successRate;
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

  void acquireWriteLock() {
    bool isLockSuccessful = pthread_rwlock_wrlock(&lock) == 0;
    if (!isLockSuccessful) {
      throw std::runtime_error("Could not acquire write lock on Memory Conscious Sampling state");
    }
  }

  void acquireReadLock() {
    bool isLockSuccessful = pthread_rwlock_rdlock(&lock) == 0;
    if (!isLockSuccessful) {
      throw std::runtime_error("Could not acquire read lock on Memory Conscious Sampling state");
    }
  }

  void releaseLock() {
    bool isUnlockSuccessful = pthread_rwlock_unlock(&lock) == 0;
    if (!isUnlockSuccessful) {
      throw std::runtime_error("Could not release lock on Memory Conscious Sampling state");
    }
  }

  void resampleWithGlobalSamplingRate() {
    for (auto it = frequency_map.begin(); it != frequency_map.end(); it++) {
      auto key = it->first;
      auto originalCount = it->second;
      if (originalCount == 0) {
        continue;
      }
      int sampleSize = getNumberSampled(originalCount, global_samplingRate);
      frequency_map[key] = sampleSize;
      if (sampleSize == 0) {
        surviving_groups--;
      }
    }
  }

  void resampleMap(double samplingRate) {
    acquireWriteLock();
    if (samplingRate >= global_samplingRate) {
      releaseLock();
      return;
    }

    global_samplingRate = samplingRate;
    resampleWithGlobalSamplingRate();
    releaseLock();
  }

  bool groupWillSurvive(HashType hash) {
    return frequency_map[hash] > minimum_group_size;
  }

  bool isKeyNew(HashType hash) const {
    return frequency_map.find(hash) == frequency_map.end();
  }

  void record_group(HashType hash) {
    if (isKeyNew(hash)) {
      frequency_map[hash] = 0;
      surviving_groups++;
    }
    frequency_map[hash]++;
  }

  void queue_resample_if_necessary(double rate) {
    if (isTooMuchMemoryUsed()) {
      double newSamplingRate = rate * reductionRate;
      resampleMap(newSamplingRate);
    }
  }

 public:
    friend class <?=$className?>;

    <?=$className?>ConstantState() {
      pthread_rwlock_init(&lock, NULL);
    }

    void AddItem(HashType hash_of_group) {
      acquireReadLock();
      double sampling_rate_during_step = global_samplingRate;
      bool isItemSampled = isBernoulliSuccess(sampling_rate_during_step);
      if (!isItemSampled) {
        releaseLock();
        return;
      }
  
      record_group(hash_of_group);
      releaseLock();
      queue_resample_if_necessary(sampling_rate_during_step);
    }

    bool IsGroupSurvivor(HashType hash) {
      return !isKeyNew(hash) && frequency_map.at(hash) > 0;
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
      'map', 'cstdlib', 'string', 'iostream'];
    $user_headers = [];
    $lib_headers  = ['base\HashFct.h'];
    $libraries    = ['armadillo'];
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
  using State = <?=$state?>;
  State& state;

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
    state.AddItem(hash);
  }

  void AddState(<?=$className?>& other)  {}

  bool IsGroupSurvivor(<?=const_typed_ref_args($inputs)?>) const {
    auto hash = chainedHash(<?=args($inputs)?>);
    return state.IsGroupSurvivor(hash);
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
