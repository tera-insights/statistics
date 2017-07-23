<?
// This version of MCS only limits the number of produced groups.
function Memory_Conscious_Sampling(array $t_args, array $inputs, array $outputs)
{
    // Class name is randomly generated
    $className = generate_name("MemoryConsciousSampling");

    // Setting output types.
    $outputs = array_combine(array_keys($outputs), $inputs);

    // The dimension of the data, i.e. how many elements are in each item.
    $dimension = count($inputs);

    // Processing of template arguments
    $minimumGroupSize = $t_args['minimumGroupSize'];
    $maximumGroupsAllowed = $t_args['maximumGroupsAllowed'];
    $initialSamplingRate = $t_args['initialSamplingRate'];
    $reductionRate = $t_args['reductionRate'];

    // Array for generating inline C++ code;
    $codeArray = array_combine(array_keys($outputs), range(0, $dimension - 1));

    $sys_headers  = ['math.h', 'armadillo', 'random', 'vector', 'stdexcept',
      'map', 'cstdlib', 'string', 'stringstream'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $extra        = [];
    $properties   = ['tuples'];
    $result_type  = ['multi'];
?>

class <?=$className?>;

class <?=$className?> {
 public:
  const int minimum_group_size = <?=$minimumGroupSize?>;
  const int maximum_groups_allowed = <?=$maximumGroupsAllowed?>;
  double samplingRate = <?=$initialSamplingRate?>;
  const double reductionRate = <?=$reductionRate?>;
  const bsl::string key_separator = "|";

  // The index used to iterate during GetNextResult;
  double return_counter;

  // Maps hashed group key to the (effective) number of tuples we have seen for
  // this group. I use the word "effective" because this algorithm resamples if
  // it sees too many groups that will survive.
  std::map<std::string, int> frequency_map;

  bool groupWillSurvive(bsl::string group) {
    return frequency_map[group] > minimum_group_size;
  }

  int numberOfGroupsThatWillSurvive() {
    int survivors = 0;
    for (auto it = frequency_map.begin(); it != frequency_map.end(); it++) {
      if (groupWillSurvive(it->first())) {
        survivors++;
      }
    }
    return survivors;
  }

  bool isTooMuchMemoryUsed() {
    return numberOfGroupsThatWillSurvive() > maximum_groups_allowed;
  }

  bool isKeyNew(bsl::string key) {
    return frequency_map.find(key) == frequency_map.end();
  }

  bool shouldSampleTuple() {
    double random = rand() * 1.0 / RAND_MAX;
    return random < samplingRate;
  }

  std::string get_group_key(<?=const_typed_ref_args($inputs)?>) {
    std::stringstream group_key_stream;
<?  foreach (array_keys($inputs) as $counter => $name) { ?>
      group_key << key_separator << <?=$name?>;
<?  } ?>
    return group_key_stream.str();
  }

  // `key` as created by `get_group_key`
  std::vector<std::string> key_to_grouping_attributes(std::string key) {
    std::vector<std::string> result;
    while (i != std::string::npos) {
      int next_index = key.substr(i).find(key_separator);
      std::string token = key.substr(i, next_index);
      result.push_back(token);
      i = next_index;

    }
    return result;
  }

  void resampleMap() {
    for (auto it = frequency_map.begin(); it != frequency_map.end(); it++) {
      int newNumberSampled = 0;
      for (int trial = 0; it->second; trial++) {
        if (shouldSampleTuple()) {
          newNumberSampled++;
        }
      }
      frequency_map[it->first] = newNumberSampled;
    }
  }

 public:
  // Constructor, nothing of note.
  <?=$className?>()
      : count(0),
        return_counter(0),
        frequency_map(0) {
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    if (!shouldSampleTuple()) {
      return;
    }

    const std::string key = get_group_key($inputs);

    if (isKeyNew(key)) {
      map[key] = 0;
    }
    map[key]++;
    if (isTooMuchMemoryUsed()) {
      currentSamplingRate *= reductionRate;
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

  const std::map<std::string, int> &GetFrequencyMap() {
    return frequency_map;
  }

  // A typical result function for a GLA of type multi. The arguments are
  // references and are then changed in the body.
  bool GetNextResult(<?=typed_ref_args($outputs)?>) {
    std::vector<std::string> groups_to_return;
    for (auto it = frequency_map.begin(); frequency_map.end(); it++) {
      if (groupWillSurvive(it->first)) {
        groups_to_return.push_back(it->first);
      }
    }

    if (return_counter >= groups_to_return.size()) {
      return false;
    }

    std::string group_key = groups_to_return.at(return_counter);
    std::vector<std::string> attributes = key_to_grouping_attributes(group_key);
<?  foreach (array_keys($outputs) as $counter => $name) { ?>
      <?=$name?> = attributes[<?=$counter?>];
<?  } ?>
    return_counter++;
    return true;
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
    ];
}
?>
