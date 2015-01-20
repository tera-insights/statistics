<?
function Big_Matrix(array $t_args, array $inputs, array $outputs)
{
    // Class name randomly generated.
    $className = generate_name("BigM");

    // Initializiation of argument names.
    $key = array_keys($inputs)[0];  // Name of the input representing key.
    $vec = array_keys($inputs)[1];  // Name of the input representing values.

    $x   = array_keys($outputs)[0];
    $y   = array_keys($outputs)[1];
    $val = array_keys($outputs)[2];

    // Setting output types.
    array_set_index($outputs, 0, $inputs[$key]);
    array_set_index($outputs, 1, $inputs[$key]);
    array_set_index($outputs, 2, lookupType("base::double"));

    // Initialization of local variables from template arguments
    $block = get_default($t_args, 'block', 40);
    $scale = get_default($t_args, 'scale', 2);
    $width = get_default($t_args, 'length', 100);
    $height = $inputs[$vec]->get('size');

    $sys_headers = ['armadillo', 'limits'];
    $user_headers = [];
    $lib_headers = [];
    $libraries = ['armadillo'];
    $extra = [];
    $result_type = 'fragment'
?>

using namespace arma;
using namespace std;

class <?=$className?>;

class <?=$className?> {
 public:
  // The side length of each blocking square.
  static const constexpr unsigned int kBlock = <?=$block?>;

  // The length of each column in the data matrix.
  static const constexpr unsigned int kHeight = <?=$height?>;

  // The initial width of the data matrix.
  static const constexpr unsigned int kWidth = <?=$width?>;

  // The proportion at which the dynamic matrix grows.
  static const constexpr unsigned int kScale = <?=$scale?>;

  struct Iterator {
    int col, num_cols, col_shift, row, num_rows, row_shift;
    bool diagonal;
    mat block;
  };

 private:
  // The data matrix being constructed item by item. The width of this matrix
  // is increased when necessary as per a dynamic array.
  mat data;

  // The set of keys processed whose indices correspond to the rows in data.
  uvec keys;

  // The variance and the mean for each item.
  rowvec stddevs, means;

  // The number of rows processed by this state.
  unsigned int count;

 public:
  <?=$className?>()
      : data(kHeight, kWidth),
        keys(kWidth),
        count(0) {
  }

  // Basic dynamic array allocation.
  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    if (count == data.n_cols) {
      data.resize(kHeight, kScale * data.n_cols);
      keys.resize(kScale * keys.n_elem);
    }
    data.col(count) = vec(<?=$vec?>.data(), kHeight);
    keys(count) = <?=$key?>;
    count++;
  }

  // Empty rows are stripped such that white space will only ever be at the end
  // of both keys and data.
  void AddState(<?=$className?> &other) {
    data.resize(kHeight, count);
    keys.resize(count);
    data = join_rows(data, other.data);
    keys = join_rows(keys, other.keys);
    count += other.count;
  }

  int GetNumFragments() {
    // The remaining whitespace is stripped.
    data.resize(kHeight, count);
    keys.resize(count);

    // Computing the relevant statistics for each item.
    stddevs = stddev(data, 1);
    means = mean(data);

    // The number of blocks is computed. For a full description of the blocking
    // scheme, refer to the top of the page.
    int ratio = (count - 1) / kBlock + 1;
    int num_fragments = ratio * (ratio + 1) / 2;
    return num_fragments;
  }

  Iterator* Finalize(int fragment) {
    // The co-ordinates of the current block in the blocking grid are computed.
    int block_col = (sqrt(1 + 8 * fragment) - 1) / 2;
    int block_row = fragment - block_col * (block_col + 1) / 2;

    // The span of this block with regard to the covariance matrix.
    int first_col = kBlock * block_col;
    int first_row = kBlock * block_row;
    int final_col = std::min(count, first_col + kBlock);
    int final_row = std::min(count, first_row + kBlock);

    // The block, i.e. the covariance submatrix, is computed.
    mat col_items = data.cols(first_col, final_col - 1);
    mat row_items = data.cols(first_row, final_row - 1).t();

    rowvec col_means =   means.subvec(first_col, final_col - 1);
    colvec row_means =   means.subvec(first_row, final_row - 1).t();
    rowvec col_stds  = stddevs.subvec(first_col, final_col - 1);
    colvec row_stds  = stddevs.subvec(first_row, final_row - 1).t();

    mat block = (row_items * col_items / kHeight - row_means * col_means)
              / (row_stds * col_stds);

    // If the block lies on the main diagonal, not all of its elements are used.
    // Because the co-variance matrix is symmetric, only the upper triangular
    // portion of the matrix is returned.
    bool diagonal = (block_col == block_row);
    return new Iterator {0, (int) col_means.n_elem, first_col,
                         0, (int) row_means.n_elem, first_row, diagonal, block};
  }

  bool GetNextResult(Iterator* it, <?=typed_ref_args($outputs)?>) {
    if (it->col == it->num_cols)
      return false;
    <?=$x?> = <?=$outputs[$x]?>(keys(it->row + it->row_shift));
    <?=$y?> = <?=$outputs[$y]?>(keys(it->col + it->col_shift));
    <?=$val?> = it->block(it->row, it->col);
    if ((it->diagonal && it->row == it->col) || it->row == it->num_rows - 1) {
      it->row = 0;
      it->col++;
    } else {
      it->row++;
    }
    return true;
  }
};

typedef <?=$className?>::Iterator <?=$className?>_Iterator;

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
