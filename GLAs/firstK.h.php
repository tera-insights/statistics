<?
// http://www.cs.umd.edu/~samir/498/vitter.pdf
function Simple_Sampling(array $t_args, array $inputs, array $outputs)
{
    // Class name is randomly generated
    $className = generate_name("RsvSamp");

    // Setting output types.
    $outputs = array_combine(array_keys($outputs), $inputs);

    // The dimension of the data, i.e. how many elements are in each item.
    $dimension = count($inputs);

    // Initialization of local variables from template arguments
    $sampleSize = $t_args['size'];

    // Array for generating inline C++ code;
    $codeArray = array_combine(array_keys($outputs), range(0, $dimension - 1));

    $sys_headers  = ['math.h', 'armadillo', 'random', 'vector', 'stdexcept'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $extra        = [];
    $properties   = ['tuples'];
    $result_type  = ['multi'];
?>

using namespace std;

class <?=$className?>;

class <?=$className?> {
 public:
  using Tuple = std::tuple<<?=typed($inputs)?>>;

  // The total amount of entries to be returned by the sampling.
  const double kSampleSize = <?=$sampleSize?>;

  // An STL vector containing the tuples of data that represent the sample. Due
  // to the small size of the sample, this can be directly stored in memory.
  vector<Tuple> sample;

  // The tuple that contains the current data whose elements are manually set.
  // This is a member variable as opposed to a local one to avoid allocating
  // memory repeatedly.
  Tuple item;

  // The total number of items processed so far.
  int count;

  // The index used to iterate during GetNextResult;
  double return_counter;

 public:
  // Constructor, nothing of note.
  <?=$className?>()
      : sample(kSampleSize),
        count(0),
        return_counter(0) {
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    if (count < kSampleSize) {
<?  foreach (array_keys($inputs) as $counter => $name) { ?>
      get<<?=$counter?>>(item) = <?=$name?>;
<?  } ?>
      sample[count] = item;
      count++;
    }
  }

  // This method is intentionally empty, as nothing has to be done to maintain
  // the sample.
  void AddState(<?=$className?>& other)  {
  }

  // Finalize is a dummy function for this GLA and only exists because it is
  // required for a GLA of type multi. Setting return_counter is only a safety
  // check as it should already be 0.
  void Finalize() {
    return_counter = 0;
  }

  // A typical result function for a GLA of type multi. The arguments are
  // references and are then changed in the body.
  bool GetNextResult(<?=typed_ref_args($outputs)?>) {
    if (return_counter >= min(kSampleSize, count))
      return false;
<?  foreach (array_keys($outputs) as $counter => $name) { ?>
    <?=$name?> = get<<?=$counter?>>(sample[return_counter]);
<?  } ?>
    return_counter++;
    return true;
  }

  // This function is intended to be called when passing this GLA as a state.
  Tuple GetSample(long counter) {
    if (counter < kSampleSize)
      return sample[counter];
    else
      throw std::invalid_argument("Sample size exceeded number of rows.");
  }

  inline const vector<Tuple>& GetTuples() {
    return sample;
  }

  int GetSize() {
    if (count < kSampleSize)
      return count;
    else
      return kSampleSize;
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
