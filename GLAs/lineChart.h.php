<?
function Line_Chart(array $t_args, array $inputs, array $outputs)
{
    // Class name randomly generated.
    $className = generate_name('LineC');

    // Initialization of local variables from template arguments.
    $length = $t_args['length']; // Number of points in the graph.
    $normalize = get_default($t_args, 'normalize', false); // Should chart be normalized.
    $p = get_default($t_args, 'p', 2);

    // Initialization of local variables from input names.
    $index = array_keys($inputs)[0];  // Index of the corresponding point.
    $value = array_keys($inputs)[1];  // Value of the corresponding point.

    $chart = array_keys($outputs)[0];

    // Setting output type.
    $output = lookupType('FixedArray',
                         ['type' => $inputs[$value], 'size' => $length]);
    array_set_index($outputs, 0, $output);

    grokit_assert(canConvert($inputs[$index], lookupType('base::integer')),
                  "Line Chart: input [0] must be convertible to integer.\n" .
                  'Received: ' . $inputs[$index]);

    grokit_assert(canConvert($inputs[$value], lookupType('base::double')),
                  "Line Chart: input [1] must be convertible to double.\n" .
                  'Received: ' . $inputs[$value]);

    $sys_headers = ['armadillo', 'limits'];
    $user_headers = [];
    $lib_headers = [];
    $libraries = ['armadillo'];
?>

using namespace arma;
using namespace std;

class <?=$className?>;

class <?=$className?> {
 public:
  // The number of points in the line chart.
  static const constexpr unsigned int kLength = <?=$length?>;

  // Epsilon used to track which indices have been filled.
  static const constexpr float kEpsilon = numeric_limits<float>::epsilon();

  // The type of the data being processed.
  using Type = <?=$inputs[$value]?>;

 private:
  // The chart being constructed, initially filled with -kEpsilon.
  Col<Type>::fixed<kLength> chart;

 public:
  <?=$className?>() {
    chart.fill(-kEpsilon);
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    if (0 <= <?=$index?> && <?=$index?> < kLength)
      chart(<?=$index?>) += <?=$value?> + kEpsilon;
  }

  void AddState(<?=$className?> &other) {
    chart += other.chart + kEpsilon;
  }

  void GetResult(<?=typed_ref_args($outputs)?>) {
<?  if ($normalize) { ?>
    chart = normalise(chart, <?=$p?>);
<?  } ?>
    <?=$chart?>.from_memory(chart.memptr());
  }

  inline Col<Type>::fixed<kLength> GetChart() const {
    return chart;
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
        'iterable'       => false,
        'input'          => $inputs,
        'output'         => $outputs,
        'result_type'    => 'single',
    ];
}
?>
