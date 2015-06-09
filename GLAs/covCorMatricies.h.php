<?
function CovCor_Matrix(array $t_args, array $inputs, array $outputs)
{
    // Setting output type
    array_set_index($outputs, 0, lookupType("base::JSON"));

    // Class name is randomly generated
    $className = generate_name("Cov");

    // The dimension of the data, i.e. how many elements are in each item.
    $dim = count($inputs);

    // Whether to produce a correlation matrix, a correlation matrix, or both.
    $which = $t_args['which'];
    $correlation = $which == "cor" || $which == "both";
    $covariance  = $which == "cov" || $which == "both";
    grokit_assert($correlation || $covariance,
        "Incorrect specification for which matrix to produce.");

    // Array for generating inline C++ code;
    $codeArray = array_combine(array_keys($inputs), range(0, $dim - 1));
?>

using namespace std;
using namespace arma;


class <?=$className?>;

class <?=$className?> {
 public:
  const int dim = <?=$dim?>;

 private:
  // A vector whose elements will be individually set to contain the current
  // item's data. The data is packaged as a vector due to the reliance on
  // armadillo. This exists as a member variable as opposed to a local variable
  // in order to avoid repeated memory allocation.
  vec::fixed<<?=$dim?>> x;

  // This is the sum of the items, used to compute the mean.
  vec::fixed<<?=$dim?>> sum;

  // The sum of x * x^T.
  mat::fixed<<?=$dim?>, <?=$dim?>> XXT;

  // The number of elements processed, used to compute the average.
  long count;

 public:
  <?=$className?>()
      : x(),
        sum(zeros<vec>(<?=$dim?>)),
        XXT(zeros<mat>(<?=$dim?>, <?=$dim?>)),
        count(0) {
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
<?  foreach ($codeArray as $name => $counter) { ?>
    x[<?=$counter?>] = <?=$name?>;
<?  } ?>
    sum += x;
    XXT += x * x.t();
    count++;
  }

  void AddState(<?=$className?> &other) {
    sum += other.sum;
    XXT += other.XXT;
    count += other.count;
  }

  void GetResult(<?=typed_ref_args($outputs)?>) {
    Json::Value _result(Json::objectValue);
<?  $vars = [ 'count' => 'count', 'mean' => 'avg', 'variance' => 'var' ];
    if ($covariance)
        $vars['covariance'] = 'cov';
    if ($correlation)
        $vars['correlation'] = 'cor';
?>

    vec::fixed<<?=$dim?>> avg(sum / count);

    vec::fixed<<?=$dim?>> var(XXT.diag() / count);
    var -= avg % avg;
    var *= 1.0 + 1.0 / (count - 1);

    mat::fixed<<?=$dim?>, <?=$dim?>> cov(XXT);
    cov -= avg * sum.t();
    cov /= count - 1;

<?  if ($correlation) { ?>
    mat::fixed<<?=$dim?>, <?=$dim?>> cor(cov);
    cor /= sqrt(var) * sqrt(var.t());
<?  } ?>

<?  ProduceResult(array_keys($outputs)[0], $vars); ?>
    cout << <?=array_keys($outputs)[0]?>.get() << endl;
  }
};

<?
    return array(
        'kind'           => 'GLA',
        'name'           => $className,
        'system_headers' => array('armadillo', 'iostream'),
	'lib_headers'    => array('ArmaJson'),
        'input'          => $inputs,
        'output'         => $outputs,
        'result_type'    => 'single',
    );
}
?>
