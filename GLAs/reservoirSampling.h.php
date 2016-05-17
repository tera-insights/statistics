<?
// http://www.cs.umd.edu/~samir/498/vitter.pdf
function Reservoir_Sampling(array $t_args, array $inputs, array $outputs)
{
    // Class name is randomly generated
    $className = generate_name("RsvSamp");

    // Setting output types.
    $outputs = array_combine(array_keys($outputs), $inputs);

    // Processing of template arguments
    $coefficient = get_default($t_args, 'coefficient', 22);
    $sampleSize = $t_args['size'];

    $sys_headers  = ['math.h', 'armadillo', 'random', 'vector', 'stdexcept'
                     'algorithm'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $extra        = [];
    $properties   = ['tuples'];
    $result_type  = ['multi'];
?>

class <?=$className?> {
 public:
  using Tuple = std::tuple<<?=typed($inputs)?>>;

  // The total amount of entries to be returned by the sampling.
  const double kSize = <?=$sampleSize?>;

  // The coefficient for the threshold at which we should switch algorithms to
  // decide skip amount. The threshold is this coefficient multiplied by the
  // sample size. Vitter showed empirically that this should be 22; this may be
  // subject to change.
  const double kThresholdCoefficient = <?=$coefficient?>;

  // The value at which we should switch algorithms to decide skip amount.
  const double kThreshold = kSize * kThresholdCoefficient;

 private:
  // An STL vector containing the tuples of data that represent the sample. Due
  // to the small size of the sample, this can be directly stored in memory.
  std::vector<Tuple> sample;

  // The tuple that contains the current data whose elements are manually set.
  // This is a member variable as opposed to a local one to avoid allocating
  // memory repeatedly.
  Tuple item;

  // The total number of items processed so far. It is a double to avoid integer
  // division in calculations.
  double count;

  // The counter for GetNextResult; denotes how many items have been returned
  // thus far. This is double because it is compared to count and casting is
  // possibly costly and should be avoided regardless of cost.
  double return_counter;

  // The total number of items processed after the reservoir is filled. Always
  // equal to count - sample size. Currently unused due to local variable being
  // instantiated in all relevant function calls.
  double number_outside_reservoir;

  // Various variables used to generate toSkip.
  double V, U, R, lhs, rhs, y, quot;

  // A counter variable representing how many more items we should skip.
  // Modified directly in the toSkipSimple functions and the toSkipGlobal
  // functions. Initialized with -1 as a safety.
  double P;

  // A standard uniform distribution which generates a double between 0
  // (inclusive) and 1 (exclusive). By Vitter's algorithm, 0 should not be
  // generated (as it would propagate an infinite loop). Hence, a sample random
  // variate is generated and subtracted from 1. Vitter's algorithm specifies
  // that 1 should also not be generated, but this is much less a problem than
  // generating 0 and is ignored.
  std::uniform_real_distribution<double> standard_uniform_distribution;

  // A uniform discrete distribution that generates a random variate between 0
  // and sample size - 1, inclusive. This is used to randomly decide which
  // element of sample should be replaced.
  std::uniform_int_distribution<long> index_distribution;

  // A random engine used to generated random variates for both of the above
  // distributions. Parallelization is not a problem for randomness as each
  // state has its own member generators.
  std::default_random_engine generator;

  // The only random variate whose value must be preserved across function
  // calls. By doing so, the number of calls to a random number generator is
  // reduced by roughly a fourth.
  double W;

 public:
  // Constructor, nothing of note.
  <?=$className?>()
      : sample(kSize),
        count(0),
        return_counter(0),
        P(-1),
        standard_uniform_distribution(0.0, 1.0),
        index_distribution(0, kSize - 1),
        generator(random_device()()),
        W(exp(log(1 - standard_uniform_distribution(generator)) / kSize)) {}

  // Simple function to insert a chosen tuple randomly into the resovoir.
  void AddCurrentItem(const Tuple& x) {
    long random_index = index_distribution(generator);
    sample[random_index] = x;
  }

  // For full explanation of the following, see Vitter's paper on reservoir
  // sampling. Modifies P via the simple algorithm described by Vitter. Local
  // variables are used to refrain from modifying count.  Modifying count
  // directly here and not when skipping elements requires less increments, but
  // not to a noticeable degree.
  void CalculateToSkipSimple() {
    double dummy_count = count + 1;
    number_outside_reservoir = dummy_count - kSize;
    V = 1.0 - standard_uniform_distribution(generator);
    quot = (double) number_outside_reservoir / (double) dummy_count;
    P = 0;  // toSkip should already be 0, merely a safety check.
    while (quot > V) {
      P++;
      dummy_count++;
      number_outside_reservoir++;
      // Quot is reduced appropiately
      quot *= number_outside_reservoir / dummy_count;
    }
  }

  // Modifies P via the simple algorithm described by Vitter. Local variables
  // are used to refrain from modifying count. As above, instructions could be
  // reduced by using add assignment for count with P. This is disregarded for
  // readability, as the reduction in instructions is neglible.
  void CalculateToSkipComplex() {
    R  = count * (W - 1.0);
    do {
      double term = count - kSize + 1;
      U = 1.0 - standard_uniform_distribution(generator);
      P = floor(R);
      // The left hand side of the inequality described by Vitter
      lhs = (U * std::pow(((count + 1) /  term), 2) * (term + P)) / (count + R);
      lhs = exp(log(lhs) / kSize);
      // The right hand side of the inequality described by Vitter
      rhs = (((count + R) / (term + P)) * term) / count;
      if (lhs < rhs) {
        W = rhs / lhs;
        break;
      }
      y  = (((U * (count + 1)) / term) * (count + P + 1)) / (count + U);
      // used for calculating the falling factorial in a loop
      double denominator, numerator_limit;
      if (kSize < P) {
        denominator = count;
        numerator_limit = term + P;
      } else {
        denominator = count - kSize + P;
        numerator_limit = count + 1;
      }

      for (double numerator = count + P; numerator >= numerator_limit;
          numerator--) {
        y *= numerator / denominator;
        denominator --;
      }
      W = exp(-log(standard_uniform_distribution(generator)) / kSize);
    } while (exp(log(y) / kSize) > (count + R) / count);
  }

  // Adds an item in the typical method of a GLA. Due to the structure of a GLA,
  // the structure of Vitter's algorithm is slightly changed. The only notable
  // difference is that skipping over data entries still requires the machine to
  // read them in; this does not result in the reduction of IO time described by
  // Vitter, but this is not an issue as IO time is extremlely fast in Grokit.
  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    if (count < kSize) {
      // State: initialization of the reservoir.
<?  foreach (array_keys($inputs) as $counter => $name) { ?>
      get<<?=$counter?>>(item) = <?=$name?>;
<?  } ?>
      sample[count] = item;
      if (count == kSize - 1)
        CalculateToSkipSimple();
    } else if (P > 0) {
      // State: skip.
      P--;
    } else {
      // State : inserting item into the reservoir and recalculating P.
<?  foreach (array_keys($inputs) as $counter => $name) { ?>
      get<<?=$counter?>>(item) = <?=$name?>;
<?  } ?>
      AddCurrentItem(item);
      if (count < kThreshold) {
        // Due to low amount of items, a simple algorithm is faster.
        CalculateToSkipSimple();
      }
      else {
        // Otherwise, it is more efficient to use the complex algorithm.
        CalculateToSkipComplex();
      }
    }
    count++;
  }

  // Merges two uniform random samples. This method is not described by Vitter,
  // who instead described a recursive call for large data sets. Instead, the
  // sample of the other state is sampled and replaces the sample of the current
  // state. This is mathematically correct due to the uniform nature of both
  // samples and the exact method involved.
  void AddState(<?=$className?>& other) {
    // A counter for which member of the current state's sample is replaced
    long replace_counter = 0;
    // The probability that a member of the other's sample is kept.
    double probability = other.count / (count + other.count);
    // Iterate over the other's sample
    for (long counter = 0; counter < kSize; counter ++) {
      double random_variable = standard_uniform_distribution(generator);
      if (random_variable <= probability)
        // The corresponding element of other's sample is kept.
        sample[replace_counter++] = other.sample[counter];
    }
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
    if (return_counter >= min(kSize, count))
      return false;
<?  foreach (array_keys($outputs) as $counter => $name) { ?>
    <?=$name?> = std::get<<?=$counter?>>(sample[return_counter]);
<?  } ?>
    return_counter++;
    return true;
  }

  // This function is intended to be called when passing this GLA as a state.
  Tuple& GetSample(long counter) {
    if (counter < kSize)
      return sample[counter];
    else
      throw std::invalid_argument("Sample size exceeded number of rows.");
  }

  inline const std::vector<Tuple>& GetTuples() {
    return sample;
  }

  std::size_t GetSize() {
    return std::min(count, kSize);
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
