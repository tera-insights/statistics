<?
function Collect(array $t_args, array $inputs, array $outputs)
{
    // Class name randomly generated.
    $className = generate_name('Collect');

    // Initialization of local variables from input names.
    $inputs_ = array_combine(['index', 'vector'], $inputs);

    grokit_assert($inputs_['index']->is('categorical'),
                  "Collect: input [0] must be categorical.\n");

    grokit_assert($inputs_['vector']->is('vector'),
                  "Collect: input [1] must be a vector.\n");

    $height = $inputs_['vector']->get('size');
    $length = $inputs_['index']->get('cardinality');
    $type = $inputs_['vector']->get('type');

    // Creating the armadillo type
    $arma = $type->is('real')
        ? $type->get('size.bytes') == 4
            ? 'fvec'
            : 'vec'
        : 'svec';


    if (count($outputs) > 0) {
      // Setting output type.
      $output = lookupType(
          'statistics::Matrix',
          ['type' => $type, 'nrow' => $height, 'ncol' => $length]
      );
      array_set_index($outputs, 0, $output);
      $outputs_ = array_combine(['result'], $outputs);
    }

    $sys_headers  = ['armadillo', 'limits'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $properties   = [];
    $extra        = ['nrow' => $height, 'ncol' => $length, 'type' => $type];
?>

using namespace arma;
using namespace std;

class <?=$className?>;

class <?=$className?> {
 public:
  // The number of columns in the matrix, equal to the cardinality of the index.
  static const constexpr unsigned int kLength = <?=$length?>;

  // The number of rows in the matrix, equal to the number of additional inputs.
  static const constexpr unsigned int kHeight = <?=$height?>;

 private:
  // The matrix to be filled.
  Mat<<?=$type?>>::fixed<kHeight, kLength> data;

 public:
  <?=$className?>() {
    data.fill(0);
  }

  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    <? if ($inputs_['vector']->get('direction') == 'row') { ?>
    data.col(index) += vector.t();
    <? } else { ?>
    data.col(index) += vector;
    <? } ?>
  }

  void AddState(<?=$className?> &other) {
    data += other.data;
  }

<?  if (count($outputs) > 0) { ?>
  void GetResult(<?=typed_ref_args($outputs_)?>) {
    result = data;
  }
<?  } ?>

  inline Mat<<?=$type?>>::fixed<kHeight, kLength> GetMatrix() const {
    return data;
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
