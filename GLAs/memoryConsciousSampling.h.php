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
  std::default_random_engine generator;

  // Maps hashed group key to the (effective) number of tuples we have seen for
  // this group. I use the word "effective" because this algorithm resamples if
  // it sees too many groups that will survive.
  Map frequency_map;

  boost::shared_mutex mutex;

  bool isTooMuchMemoryUsed() {
    return frequency_map.size() > maximum_groups_allowed;
  }

  void resampleWithReductionRate() {
    auto it = frequency_map.begin();
    while (it != frequency_map.end()) {
      auto key = it->first;
      auto originalCount = it->second;
      if (originalCount == 0) {
        continue;
      }
      std::binomial_distribution<int> distribution(originalCount,
        reductionRate);
      int newSampleSize = distribution(generator);
      if (newSampleSize == 0) {
        it = frequency_map.erase(it);
      } else {
        frequency_map[key] = newSampleSize;
      }
      if (it != frequency_map.end()) {
        it++;
      }
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
    while (isTooMuchMemoryUsed()) {
      double newSamplingRate = rate * reductionRate;
      resampleMap(newSamplingRate);
      std::cout << "map size = " << frequency_map.size() << std::endl;
    }
  }

 public:
    friend class <?=$className?>;

    <?=$className?>ConstantState()
      : global_samplingRate(<?=$initialSamplingRate?>),
        frequency_map{},
        generator(RandInt()) {
    }

    void AddMap(const Map &other_map) {
      boost::unique_lock<boost::shared_mutex> lock(mutex);
      double sampling_rate_during_step = global_samplingRate;
      for (auto it = other_map.begin(); it != other_map.end(); it++) {
        std::binomial_distribution<int> distribution(it->second,
          sampling_rate_during_step);
        int sampleSize = distribution(generator);
        if (sampleSize > 0) {
          frequency_map[it->first] += sampleSize;
        }
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

  void AddState(const <?=$className?>& other) {
    auto other_map
    for (auto it = other_map.begin(); it != other_map.end(); it++) {
      frequency_map[it->first] += it->second;
    }
  }

  bool IsGroupSurvivor(<?=const_typed_ref_args($inputs)?>) const {
    auto hash = chainedHash(<?=args($inputs)?>);
    return state.IsGroupSurvivor(hash);
  }

  void ChunkBoundary() {
    state.AddMap(frequency_map);
    frequency_map.clear();
  }

  Map GetMap() const {
    return frequency_map;
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
        'chunk_boundary' => true,
        'generated_state' => $state,
    ];
}
?>
