<?
function Gather(array $t_args, array $inputs, array $outputs)
{
    // Class name is randomly generated.
    $className = generate_name('Gather');

    // Initialization of local variables from template arguments.
    $scale = get_default($t_args, 'scale',  2);
    $width = get_default($t_args, 'length', 100);

    // Initialization of local variables from inputs.
    $inputs_ = array_combine(['vector', 'tuple'], $inputs);

    grokit_assert($inputs_['vector']->is('vector'),
                  "Collect: input [0] must be a vector.\n");

    $height = $inputs_['vector']->get('size');
    $type   = $inputs_['vector']->get('type');
    $tuple  = $inputs_['tuple'];

    $sys_headers  = ['armadillo'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $properties   = ['matrix', 'tuples'];
    $extra        = [];
?>

using namespace arma;
using namespace std;

class <?=$className?>;

class <?=$className?> {
 public:
  // The type of the item vectors.
  using Inner = <?=$type?>;

  // The type of the extra tuples.
  using Tuple = <?=$tuple?>;

  // The length of each column in the data matrix.
  static const constexpr unsigned int kHeight = <?=$height?>;

  // The initial width of the data matrix.
  static const constexpr int kWidth = <?=$width?>;

  // The proportion at which the dynamic matrix grows.
  static const constexpr double kScale = <?=$scale?>;

 private:
  // The data matrix being constructed item by item. The width of this matrix
  // is increased when necessary as per a dynamic array.
  Mat<Inner> items;

  // Extra attributes that either cannot fit into the matrix or are not needed
  // in such a format.
  vector<Tuple> extra;

  // The number of inputs processed by this state.
  long count;

 public:
  <?=$className?>()
      : items(kHeight, kWidth),
        extra(),
        count(0) {
  }

  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    if (count == items.n_cols)
      items.resize(kHeight, kScale * items.n_cols);
    items.col(count) = vector;
    extra.push_back(tuple);
    count++;
  }

  void AddState(<?=$className?> &other) {
    items.resize(kHeight, count);
    items.insert_cols(count, other.items);
    extra.reserve(extra.size() + other.extra.size());
    extra.insert(extra.end(), other.extra.begin(), other.extra.end());
    count += other.count;
  }

  void FinalizeState() {
    // The remaining whitespace is stripped.
    items.resize(kHeight, count);
    cout << "<?=$className?> processed " << count << " tuples" << endl;
  }

  const Mat<Inner>& GetMatrix() const {
    return items;
  }

  const vector<Tuple>& GetTuples() const {
    return extra;
  }

  long GetCount() const {
    return count;
  }
};

<?
    return [
        'kind'              => 'GLA',
        'name'              => $className,
        'system_headers'    => $sys_headers,
        'user_headers'      => $user_headers,
        'lib_headers'       => $lib_headers,
        'libraries'         => $libraries,
        'iterable'          => false,
        'input'             => $inputs,
        'output'            => $outputs,
        'result_type'       => 'single',
        'finalize_as_state' => true,
        'properties'        => $properties,
        'extra'             => $extra,
    ];
}
?>
