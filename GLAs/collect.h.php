/**
 * This GLA takes in tuples and combines them into a matrix. 
 * $inputs must contain a pair of labels for the index and vector (in that 
 * order). 
 * $outputs contains the types of the output. 
 * 
 * As an example, consider the CSV file whose contents are:
 * 
 * 0,1 2
 * 1,3 4
 * 3,5 6
 *
 * The following code creates a matrix with values:
 *
 * 1 2
 * 3 4
 * 5 6
 *
 * data <- ReadCSV("filename.csv", c(Row = base::int, 
 *   Data = statistics::vector(size = 2, type = base::int)), header = FALSE)
 * state <- Collect(data, c(Row, Data), Matrix, size = 3)
 */

<?
function Collect(array $t_args, array $inputs, array $outputs)
{
    // Class name is randomly generated.
    $className = generate_name('Collect');

    // Initialization of local variables from input names.
    $inputs_ = array_combine(['index', 'vector'], $inputs);

    grokit_assert(   $inputs_['index']->is('categorical')
                  || array_key_exists('size', $t_args),
                  "Collect: input [0] must be categorical.\n");

    grokit_assert($inputs_['vector']->is('vector'),
                  "Collect: input [1] must be a vector.\n");

    $height = $inputs_['vector']->get('size');
    $length = $inputs_['index']->is('categorical')
        ? $inputs_['index']->get('cardinality')
        : $t_args['size'];

    $type = $inputs_['vector']->get('type');
    if ($type->is('categorical'))
        $type = lookupType('int');

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
    $properties   = ['matrix'];
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

  // The type of the matrix being constructed.
  using Matrix = Mat<<?=$type?>>::fixed<kHeight, kLength>;

 private:
  // The matrix to be filled.
  Matrix data;

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

  inline const Matrix& GetMatrix() const {
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
