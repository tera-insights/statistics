
<?
function Memory_Conscious_Hashing(array $t_args, array $inputs, array $outputs)
{
    // Randomly generate class name
    $className = generate_name("MemoryConsciousHashing");

    // Set output types
    $outputs = array_combine(array_keys($outputs), $inputs);

    // Process template arguments
    $initialThreshold = $t_args['initialThreshold'];
    $thresholdGrowthFactor = $t_args['thresholdGrowthFactor'];
    $maxProducedGroups = $t_args['maxProducedGroups'];
    $dieProbability = $t_args['dieProbability'];

    $isDebugMode = get_default($t_args, 'debug', 0) > 0;
    $tuplesPerFragment = 200000;

    $sys_headers  = ['math.h', 'armadillo', 'random', 'vector', 'stdexcept',
      'map', 'cstdlib', 'string', 'iostream'];
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
  using Map = std::map<int, int>;

  // For the fragmented result
  using Iterator = struct {
    Map::const_iterator current, end;
    std::size_t index;
  };

  const double initialThreshold = <?=$initialThreshold?>;
  const double thresholdGrowthFactor = <?=$thresholdGrowthFactor?>;
  const int maxProducedGroups = <?=$maxProducedGroups?>;
  const double dieProbability = <?=$dieProbability?>;
  static constexpr std::size_t kFragmentSize = <?=$tuplesPerFragment?>;

  // Maps hashed group key to the bucket score. Multiple groups hash to the
  // same bucket.
  Map bucket_scores;

  // Fields used for the fragment result type.
  std::vector<Map::const_iterator> iterators;

  int get_bucket_key(KeySet group) {
    
  }

  bool groupWillSurvive(KeySet group) {
    return frequency_map[group] > minimum_group_size;
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
    return numberOfGroupsThatWillSurvive() > maximum_groups_allowed;
  }

  bool isKeyNew(KeySet key) {
    return frequency_map.find(key) == frequency_map.end();
  }

  bool shouldSampleTuple() {
    double random = rand() * 1.0 / RAND_MAX;
    return random < samplingRate;
  }

  KeySet get_group_key(<?=const_typed_ref_args($inputs)?>) {
    return KeySet(<?=args($inputs)?>);
  }

  void resampleMap() {
    for (auto it = frequency_map.begin(); it != frequency_map.end(); it++) {
      auto key = it->first;
      auto originalCount = it->second;
      int newNumberSampled = 0;
      for (int trial = 0; trial < originalCount; trial++) {
        if (shouldSampleTuple()) {
          newNumberSampled++;
        }
      }
      if (newNumberSampled == 0) {
        frequency_map.erase(key);
      } else {
        frequency_map[key] = newNumberSampled;
      }
    }
  }

 public:
  // Constructor, nothing of note.
  <?=$className?>()
      : bucket_scores() {
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    if (!shouldSampleTuple()) {
      return;
    }

    const KeySet key = get_group_key(<?=args($inputs)?>);
    if (isKeyNew(key)) {
      frequency_map[key] = 0;
    }
    frequency_map[key]++;
    if (isTooMuchMemoryUsed()) {
      samplingRate *= reductionRate;
      resampleMap();
    }
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
