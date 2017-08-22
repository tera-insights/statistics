
<?
// This version of MCS only limits the number of produced groups.
function Memory_Conscious_Sampling(array $t_args, array $inputs, array $outputs)
{
    $className = generate_name("MemoryConsciousSampling");
    $outputs = ["hashed_key" => lookupType("BASE::INT")];

    // The dimension of the data, i.e. how many elements are in each item.
    $dimension = count($inputs);

    // Processing of template arguments
    $minimumGroupSize = $t_args['minimumGroupSize'];
    $maximumGroupsAllowed = $t_args['maximumGroupsAllowed'];
    $initialSamplingRate = $t_args['initialSamplingRate'];
    $reductionRate = $t_args['reductionRate'];

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

class <?=$className?> {
 public:
  using HashType = uint64_t;
  using Map = std::map<HashType, int>;

  const int minimum_group_size = <?=$minimumGroupSize?>;
  const int maximum_groups_allowed = <?=$maximumGroupsAllowed?>;
  double samplingRate = <?=$initialSamplingRate?>;
  const double reductionRate = <?=$reductionRate?>;

  // Maps hashed group key to the (effective) number of tuples we have seen for
  // this group. I use the word "effective" because this algorithm resamples if
  // it sees too many groups that will survive.
  Map frequency_map;

  uint64_t surviving_groups;

  bool groupWillSurvive(HashType hash) {
    return frequency_map[hash] > minimum_group_size;
  }

  int numberOfGroupsThatWillSurvive() {
    int survivors = 0;
    for (auto it = frequency_map.begin(); it != frequency_map.end(); it++) {
      if (groupWillSurvive(it->first)) {
        survivors++;
      }
    }
    return survivors;
  }

  bool isTooMuchMemoryUsed() {
    return surviving_groups > maximum_groups_allowed;
  }

  bool isKeyNew(HashType hash) const {
    return frequency_map.find(hash) == frequency_map.end();
  }

  bool shouldSampleTuple() {
    double random = rand() * 1.0 / RAND_MAX;
    return random < samplingRate;
  }

  HashType chainedHash(<?=const_typed_ref_args($inputs)?>) const {
    HashType result = 0;
<?  foreach ($inputs as $index => $value) {?>
    result = SpookyHash(Hash(<?=$index?>), result);
<?}?>
    return result;
  }

  void FinalizeState() {
    while (isTooMuchMemoryUsed()) {
      samplingRate *= reductionRate;
      resampleMap();
    }
  }

  void resampleMap() {
    for (auto it = frequency_map.begin(); it != frequency_map.end(); it++) {
      auto key = it->first;
      auto originalCount = it->second;
      if (originalCount == 0) {
        continue;
      }
      int newNumberSampled = 0;
      for (int trial = 0; trial < originalCount; trial++) {
        if (shouldSampleTuple()) {
          newNumberSampled++;
        }
      }
      frequency_map[key] = newNumberSampled;
      if (newNumberSampled == 0) {
        surviving_groups--;
      }
    }
  }

 public:
  <?=$className?>()
      : frequency_map(),
      surviving_groups(0) {
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    if (!shouldSampleTuple()) {
      return;
    }

    const HashType hash = chainedHash(<?=args($inputs)?>);
    if (isKeyNew(hash)) {
      frequency_map[hash] = 0;
      surviving_groups++;
    }
    frequency_map[hash]++;
  }

  // Merges `this->frequency_map` and `other->frequency_map`.
  void AddState(<?=$className?>& other)  {
    auto other_map = other.GetFrequencyMap();
    for (auto it = other_map.begin(); it != other_map.end(); it++) {
      auto key = it->first;
      if (isKeyNew(key)) {
        frequency_map[key] = 0;
        surviving_groups++;
      }
      frequency_map[key] += it->second;
    }
  }

  const std::map<HashType, int> &GetFrequencyMap() {
    return frequency_map;
  }

  bool IsGroupSurvivor(<?=const_typed_ref_args($inputs)?>) const {
    auto hash = chainedHash(<?=args($inputs)?>);
    return !isKeyNew(hash) && frequency_map.at(hash) > 0;
  }

  double GetSamplingRate() const {
    return samplingRate;
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
    ];
}
?>
