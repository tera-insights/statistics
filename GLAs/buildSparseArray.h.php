<?
// This GLA builds a fixed-sized multi-dimensional array. Given a list of scalar
// indices followed by a scalar or an array, it accumulates the data into an
// array of higher dimension.

// Template Args:
// Scale: The scaling factor the dynamically allocated matrix for the input.
// Width: The input matrix is initially allocated to hold this many inputs.
function Build_Sparse_Array(array $t_args, array $inputs, array $outputs)
{
    // Class name is randomly generated.
    $className = generate_name('BuildSparseArray');

    // Initializiation of argument names.
    $inputs_ = array_combine(['key', 'vector'], $inputs);
    $key = $inputs_['key'];

    // Information about the vector type.
    $size = $inputs_['vector']->get('size');
    $type = $inputs_['vector']->get('type');

    // Processing of template arguments.
    $scale = get_default($t_args, 'scale',  2);
    $width = get_default($t_args, 'length', 100);
    $type  = get_default($t_args, 'type',   $type);

    // diag is converted to avoid PHP boolean printing issues.
    $diag = intval($diag);

    grokit_assert(is_datatype($type),
                  "Build Sparse Array: 'type' ($type) is not a datatype.");

    $sys_headers  = ['armadillo', 'limits'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $extra        = [];
    $result_type  = ['state'];
?>

using namespace std;

class <?=$className?>;

class <?=$className?> {
 public:
  // The type of result.
  using Type = <?=$type?>;

  // The type of the indices.
  using Key = <?=$key?>;

  // The length of each column in the data matrix.
  static const constexpr int kHeight = <?=$size?>;

  // The initial width of the data matrix.
  static const constexpr int kWidth = <?=$width?>;

  // The proportion at which the dynamic matrix grows.
  static const constexpr int kScale = <?=$scale?>;

 private:
  // The data matrix being constructed item by item. The width of this matrix
  // is increased when necessary as per a dynamic array.
  Mat<Type> data;

  // The set of keys processed whose indices correspond to the rows in data.
  Col<Key> keys;

  // The number of rows processed by this state.
  int count;

 public:
  <?=$className?>()
      : data(kHeight, kWidth),
        keys(kWidth),
        count(0) {
  }

  // Basic dynamic array allocation.
  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    if (count == data.n_cols) {
      data.resize(kHeight, kScale * data.n_cols);
      keys.resize(kScale * keys.n_elem);
    }
    data.col(count) = Col<Type>(vector.data(), kHeight);
    keys(count) = key;
    count++;
  }

  // Empty rows are stripped such that white space will only ever be at the end
  // of both keys and data.
  void AddState(<?=$className?> &other) {
    data.resize(kHeight, count);
    keys.resize(count);
    data.insert_cols(count, other.data);
    keys.insert_rows(count, other.keys);
    count += other.count;
  }

  void FinalizeState() {
    data.resize(kHeight, count);
  }

  const Col<Key>& GetKeys() const {
    return keys;
  }

  const Mat<Type>& GetData() const {
    return data;
  }

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

function get_dimension($type) {
    return count(get_dimensions($type));
}
?>
