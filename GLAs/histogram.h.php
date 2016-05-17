<?
function Histogram(array $t_args, array $inputs, array $outputs)
{
    // Class name is randomly generated.
    $className = generate_name("Hist");

    // Processing of template arguments.
    $length = $t_args['length'];
    $normalize = get_default($t_args, 'normalize', false);
    $p = get_default($t_args, 'p', 1);

    // Initialization of local variables from input names.
    $index = array_keys($inputs)[0];  // Index of the corresponding point.

    grokit_assert(canConvert($inputs[$index], lookupType('base::integer')),
                  "Line Chart: input [0] must be convertible to integer.\n" .
                  'Received: ' . $inputs[$index]);

    // Setting output type
    $arrayInnerType = lookupType('base::Integer');
    $arrayType = lookupType('base::FixedArray',
                            ['type' => $arrayInnerType, 'size' => $length]);
    array_set_index($outputs, 0, $arrayType);

    $sys_headers = ['armadillo'];
    $user_headers = [];
    $lib_headers = [];
    $libraries = ['armadillo'];
?>

using namespace arma;
using namespace std;

class <?=$className?>;

class <?=$className?> {
 public:
  // The number of bins in the histogram.
  static const constexpr unsigned int kLength = <?=$length?>;

 private:
  // A vector containing the bin counts.
  ivec::fixed<kLength> bins;

 public:
  <?=$className?>() {
    bins.fill(0);
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    if (<?=$index?> >= 0 && <?=$index?> < kLength)
      bins(<?=$index?>)++;
  }

  void AddState(<?=$className?> &other) {
      bins += other.bins;
  }

  void GetResult(<?=typed_ref_args($outputs)?>) {
<?  if ($normalize) { ?>
    bins = normalise(bins, <?=$p?>);
<?  } ?>
    <?=array_keys($outputs)[0]?>.from_memory(bins.memptr());
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
