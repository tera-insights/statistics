<?
// This GLA is used compute pair-wise statistics between input vectors, such as
// correlation and covariance. This is done via the following:
// 1. Input vectors are collected into a matrix.
// 2. The output matrix is partitioned into blocks.
// 3. Each block is computed separately and outputs a list of paired stastics.
// The blocking is performed by partitioning the n inputs into k intervals
// within the matrix. Each block is then given two intervals, 0 <= k1 <= k2 < k,
// which represesent the two components of each pair-wise statistics.

// Template Args:
// block: The side length of a block, i.e. the length of each interval.
// scale: The scaling factor the dynamically allocated matrix for the input.
// width: The input matrix is initially allocated to hold this many inputs.
// type:  The type used to perform calculations. The input data type by default.
// diag:  Should diagonal entries be returned. Usually, these entries hold no
//   meaningful information, e.g. the correlation is always 1.
function Big_Matrix($t_args, $inputs, $outputs)
{
    // Class name is randomly generated.
    $className = generate_name('BigM');

    // Initializiation of argument names.
    $inputs_ = array_combine(['key', 'vector'], $inputs);

    // Information about the vector type.
    $size = $inputs_['vector']->get('size');
    $type = $inputs_['vector']->get('type');

    // Processing of template arguments.
    $block = get_default($t_args, 'block',  40);
    $scale = get_default($t_args, 'scale',  2);
    $width = get_default($t_args, 'length', 100);
    $type  = get_default($t_args, 'type',   $type);
    $diag  = get_default($t_args, 'diag',   true);

    // diag is converted to avoid PHP boolean printing issues.
    $diag = intval($diag);

    grokit_assert(is_datatype($type),
                  "BigMatrix: 'type' ($type) is not a datatype.");

    // Construction of outputs.
    $key = $inputs_['key'];
    $outputs_ = ['x' => $key, 'y' => $key, 'val' => lookupType('double')];
    $outputs = array_combine(array_keys($outputs), $outputs_);

    $sys_headers  = ['armadillo', 'limits'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $extra        = [];
    $result_type  = ['fragment'];
?>

using namespace arma;
using namespace std;

class <?=$className?>;

class <?=$className?> {
 public:
  // The side length of each blocking square.
  static const constexpr int kBlock = <?=$block?>;

  // The length of each column in the data matrix.
  static const constexpr int kHeight = <?=$size?>;

  // The initial width of the data matrix.
  static const constexpr int kWidth = <?=$width?>;

  // The proportion at which the dynamic matrix grows.
  static const constexpr int kScale = <?=$scale?>;

  // Whether the diagonal entries are kept.
  static const constexpr bool kDiag = <?=$diag?>;

  // The type of the data being processed.
  using Type = <?=$type?>;

  // The type of the indices.
  using Key = <?=$key?>;

  struct Iterator {
    // The indices of the big matrix where the upper left of this fragment is.
    int col_shift, row_shift;

    // The data associated with this fragment.
    mat block;

    // The number of columns and rows in this fragment.
    int n_cols, n_rows;

    // Whether the corresponding fragment is on the main diagonal.
    bool diagonal;

    // Used to iterate over this fragment during output.
    int col, row;

    Iterator(int col_shift, int row_shift, bool diagonal, mat&& block)
        : col_shift(col_shift),
          row_shift(row_shift),
          block(block),
          n_cols(block.n_cols),
          n_rows(block.n_rows),
          diagonal(diagonal),
<?  if ($diag) { ?>
          col(0),
<?  } else { ?>
          col(diagonal),
<?  } ?>
          row(0) {
    }
  };

 private:
  // The data matrix being constructed item by item. The width of this matrix
  // is increased when necessary as per a dynamic array.
  mat data;

  // The set of keys processed whose indices correspond to the rows in data.
  std::vector<Key> keys;

  // The variance and the mean for each item.
  arma::rowvec stddevs, means;

  // The number of rows processed by this state.
  long count;

 public:
  <?=$className?>()
      : data(kHeight, kWidth),
        count(0) {
  }

  // Basic dynamic array allocation.
  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    if (count == data.n_cols)
      data.resize(kHeight, kScale * count);
    data.col(count) = conv_to<colvec>::from(Col<Type>(vector.data(), kHeight));
    keys.push_back(key);
    count++;
  }

  // Empty rows are stripped such that white space will only ever be at the end
  // of both keys and data.
  void AddState(<?=$className?> &other) {
    data.resize(kHeight, count);
    data.insert_cols(count, other.data);
    keys.insert(keys.end(), other.keys.begin(), other.keys.end());
    count += other.count;
  }

  int GetNumFragments() {
    // The remaining whitespace is stripped.
    data.resize(kHeight, count);

    // Computing the relevant statistics for each item.
    stddevs = arma::stddev(data, 1);
    means = arma::mean(data);

    arma::uvec invalid = arma::find(stddevs == 0);
    std::cout << "Number of 0 stddevs: " << invalid.n_elem << std::endl;
    if (invalid.n_elem > 0) {
      invalid = invalid.subvec(0, std::min(5, (int) invalid.n_elem));
      char buffer [100];
      std::cout << "Keys for illegal vectors:" << std::endl;
      for (auto i : invalid) {
        keys[i].ToString(buffer);
        std::cout << " " << buffer;
      }
      std::cout << std::endl;
    }

    // The number of blocks is computed. For a full description of the blocking
    // scheme, refer to the top of the page.
    int ratio = (count - 1) / kBlock + 1;
    int n_fragments = ratio * (ratio + 1) / 2;
    return n_fragments;
  }

  Iterator* Finalize(int fragment) {
    // The co-ordinates of the current block in the blocking grid are computed.
    long block_col = (sqrt(1 + 8 * fragment) - 1) / 2;
    long block_row = fragment - block_col * (block_col + 1) / 2;

    // The span of this block with regard to the covariance matrix.
    // Both intervals are closed on the left and open on the right.
    long first_col = kBlock * block_col;
    long first_row = kBlock * block_row;
    long final_col = std::min(count, first_col + kBlock);
    long final_row = std::min(count, first_row + kBlock);

    // The dimension of the rectangular block.
    long num_cols = final_col - first_col;
    long num_rows = final_row - first_row;

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

    return new Iterator(first_col, first_row, diagonal, std::move(block));
  }

  bool GetNextResult(Iterator* it, <?=typed_ref_args($outputs_)?>) {
    if (it->col == it->n_cols)
      return false;
    x = keys[it->row + it->row_shift];
    y = keys[it->col + it->col_shift];
    val = it->block(it->row, it->col);
<?  if ($diag) { ?>
    if ((it->diagonal && it->row == it->col) || it->row == it->n_rows - 1) {
<?  } else { ?>
    if ((it->diagonal && it->row == it->col - 1) || it->row == it->n_rows - 1) {
<?  } ?>
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
