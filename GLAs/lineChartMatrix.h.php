<?
function Line_Chart_Matrix(array $t_args, array $inputs, array $outputs)
{
    // Class name randomly generated
    $className = generate_name("LineCM");

    // Initializiation of argument names.
    $key   = array_keys($inputs)[0];  // Name of the input representing key.
    $index = array_keys($inputs)[1];  // Name of the input representing time.
    $value = array_keys($inputs)[2];  // Name of the input representing value.

    // Initialization of local variables from template arguments
    $length = $t_args['length'];  // Number of points in the graph per key.

    grokit_assert($inputs[$key]->is("categorical"),
                  "Line Chart Matrix: input [0] must be categorical.");
    grokit_assert(canConvert($inputs[$index], lookupType('base::integer')),
                  'Line Chart Matrix: input [1] must be convertible to integer.\n' .
                  'Received: ' . $inputs[$index]);
    grokit_assert(canConvert($inputs[$value], lookupType("base::float")),
                  "Line Chart Matrix: input [2] must be convertible to double.\n" .
                  "Received: " . $inputs[$value]);

    // Setting output type
    if (count($outputs) > 0)
        array_set_index($outputs, 0, lookupType(
            'base::array',
            ['size' => $length, 'type' => lookupType("base::double")]
        ));

    // Cardinality of the key.
    $cardinality = $inputs[$key]->get("cardinality");

    $sys_headers = ['armadillo', 'limits'];
    $user_headers = [];
    $lib_headers = [];
    $libraries = ['armadillo'];
    $extra = ['key' => $inputs[$key]];
    $result_type = (count($outputs) > 0) ? 'multi' : 'state';
?>

using namespace arma;
using namespace std;

class <?=$className?>;

class <?=$className?> {
 public:
  // The number of points in each line chart.
  static const constexpr unsigned int kLength = <?=$length?>;

  // The number of charts to be constructed.
  static const constexpr unsigned int kCardinality = <?=$cardinality?>;

  // Epsilon used to track which indices have been filled.
  static const constexpr float kEpsilon = numeric_limits<float>::epsilon();

 private:
  // The charts being constructed, initially filled with -kEpsilon.
  mat::fixed<kLength, kCardinality> charts;

<?  if (count($outputs) > 0) { ?>
  // The counter used to iterate over charts in GetNextResult.
  int counter;
<?  } ?>

 public:
  <?=$className?>() {
    charts.fill(-kEpsilon);
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    if (0 <= <?=$index?> && <?=$index?> < kLength)
      charts(<?=$index?>, <?=$key?>) += <?=$value?> + kEpsilon;
  }

  void AddState(<?=$className?> &other) {
    charts += other.charts + kEpsilon;
  }

<?  if (count($outputs) > 0) { ?>
  void Finalize() {
    counter = 0;
  }

  bool GetNextResult(<?=typed_ref_args($outputs)?>) {
    if (counter == kCardinality)
      return false;
    <?=array_keys($outputs)[0]?>.from_memory(charts.colptr(counter));
    counter++;
    return true;
  }
<?  } ?>

  inline mat::fixed<kLength, kCardinality> GetCharts() const {
    return charts;
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
