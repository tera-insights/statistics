<?
function Summary(array $t_args, array $inputs, array $outputs)
{
    // Class name randomly generated
    $className = generate_name("Summ");

    // The dimension of the data, i.e. how many elements are in each item.
    $dimension = count($inputs);

    // Array for generating inline C++ code;
    $codeArray = array_combine(array_keys($inputs), range(0, $dimension - 1));

    // Setting output type
    for ($counter = 0; $counter < 5 * $dimension; $counter++)
      array_set_index($outputs, $counter, lookupType("base::DOUBLE"));
    array_set_index($outputs, $counter, lookupType("base::INT"));

    $names = array_keys($outputs);
    $variables = ["mean", "std_dev", "range", "min", "max"];
    $count = array_pop($names);

    $sys_headers = ['armadillo', 'iostream', 'limits'];
    $user_headers = [];
    $lib_headers = [];
?>

using namespace std;
using namespace arma;

class <?=$className?>;

class <?=$className?> {
 private:
  // A vector whose elements will be individually set to contain the current
  // item's data. The data is packaged as a vector due to the reliance on
  // armadillo. This exists as a member variable as opposed to a local variable
  // in order to avoid repeated memory allocation.
  vec::fixed<<?=$dimension?>> item;

  // This is the sum of the items, used to compute the means.
  vec::fixed<<?=$dimension?>> sum;

  // This is the sum of the items squared, used to compute the variances.
  vec::fixed<<?=$dimension?>> sum_squared;

  // These contain the min and max for each input, respectively.
  vec::fixed<<?=$dimension?>> min, max;

  // The number of elements processed, used to compute the means.
  long count;

 public:
  <?=$className?>()
      : item(),
        count(0),
        sum(zeros<vec>(<?=$dimension?>)),
        sum_squared(zeros<vec>(<?=$dimension?>)) {
    max.fill(-numeric_limits<double>::infinity());
    min.fill( numeric_limits<double>::infinity());
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
<?  foreach ($codeArray as $name => $counter) { ?>
    item[<?=$counter?>] = <?=$name?>;
<?  } ?>
    if (any(item != item))
      cout << "nan seen" << endl;
    sum += item;
    sum_squared += item % item;
    min = arma::min(min, item);
    max = arma::max(max, item);
    count++;
  }

  void AddState(<?=$className?> &other) {
    sum += other.sum;
    sum_squared += other.sum_squared;
    min = arma::min(min, other.min);
    max = arma::max(max, other.max);
    count += other.count;
  }

  void GetResult(<?=typed_ref_args($outputs)?>) {
    vec::fixed<<?=$dimension?>> mean(sum / count);
    vec::fixed<<?=$dimension?>> var(sum_squared / count - mean % mean);
    var *= 1.0 + 1.0 / (count - 1);
    vec::fixed<<?=$dimension?>> std_dev(sqrt(var));
    vec::fixed<<?=$dimension?>> range(max - min);
    cout << "sum" << endl << sum << endl;
    cout << "count" << endl << count << endl;
    cout << "mean" << endl << mean << endl;
<?  foreach ($names as $counter => $output) { ?>
    <?=$output?> = <?=$variables[$counter % 5]?>(<?=floor($counter / 5)?>);
<?  } ?>
    <?=$count?> = count;
  }
};

<?
    return [
        'kind'           => 'GLA',
        'name'           => $className,
        'system_headers' => $sys_headers,
        'user_headers'   => $user_headers,
        'lib_headers'    => $lib_headers,
        'input'          => $inputs,
        'output'         => $outputs,
        'result_type'    => 'single',
    ];
}
?>
