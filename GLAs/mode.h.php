<?
function Mode($t_args, $inputs, $outputs)
{
    // Class name randomly generated.
    $className = generate_name('Collect');

    // Initialization of local variables from input names.
    $inputs_ = array_combine(['index'], $inputs);

    grokit_assert($inputs_['index']->is('categorical'),
                  "Mode: input [0] must be categorical.\n");

    $length = $inputs_['index']->get('cardinality');

    $outputs_ = array_combine(['result'], $inputs);
    $outputs  = array_combine(array_keys($outputs), $outputs_);

    $sys_headers  = ['armadillo'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $properties   = [];
    $extra        = [];
?>

using namespace arma;
using namespace std;

class <?=$className?>;

class <?=$className?> {
 public:
  // The number of columns in the matrix, equal to the cardinality of the index.
  static const constexpr unsigned int kLength = <?=$length?>;

 private:
  // The vector containing the counts.
  vec::fixed<kLength> counts;

 public:
  <?=$className?>() {
    counts.fill(0);
  }

  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    counts(index)++;
  }

  void AddState(<?=$className?> &other) {
    counts += other.counts;
  }

  void GetResult(<?=typed_ref_args($outputs_)?>) {
    counts.max(result);
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
        'properties'     => $properties,
        'extra'          => $extra,
    ];
}
?>
