<?
function Big_Matrix_Batch($t_args, $outputs, $states)
{
    // Class name is randomly generated.
    $className = generate_name('BigMatrixBatch');

    $numStates = count($states);
    grokit_assert(in_array($numStates, [1, 2]),
                  "BMB: Illegal number of states: $numStates.");

    $states_ = array_values($states);
    $type = array_get_index($states_[0]->inputs(), 0)->get('type');
    $type = array_get_index($states_[0]->inputs(), 0)->get('size');
    $key =  array_get_index($states_[0]->inputs(), 1);

    $single = count($states) == 1;

    // Initialization of local variables from template arguments.
    $block = get_default($t_args, 'block',  40);
    $diag  = get_default($t_args, 'diag', true);

    // Setting output types.
    $outputs_ = [$key, $key, $type];
    $outputs = array_combine(array_keys($outputs), $outputs_);

    $sys_headers  = ['armadillo', 'vector'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $extra        = [];
    $result_type  = ['fragment'];
?>

using namespace arma;
using namespace std;

class EndGLA {
 public:
  EndGLA() {}

  void AddState(AnswerGLA other) {}

  bool ShouldIterate() {
    return false;
  }
};

class <?=$className?>;

class <?=$className?> {
 public:
  struct Iterator {
    // The indices of the big matrix where the upper left of this fragment is.
    int col_shift, row_shift;

    // The data associated with this fragment.
    Mat<Type> block;

    // The number of columns and rows in this fragment.
    int n_cols, n_rows;

    // Whether the corresponding fragment is on the main diagonal.
    bool diagonal;

    // Used to iterate over this fragment during output.
    int col, row;

    Iterator(int col_shift, int row_shift, int n_cols, int n_rows,
             bool diagonal, Mat<Type>&& block)
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

  struct Task {};

  struct LocalScheduler {
    LocalScheduler() {}

    bool GetNextTask(Task& task) {
      return false;
    }
  };

  // The type of the extra tuples.
  using Tuple = <?=$tuple?>;

  // The inner GLA being used.
  using cGLA = EndGLA;

  // The type of the workers.
  using WorkUnit = pair<LocalScheduler*, cGLA*>;

  // The type of the container for the workers.
  using WorkUnits = vector<WorkUnit>;

  // The type of result.
  using Type = <?=$type?>;

  // The type of the indices.
  using Key = <?=$key?>;

  // The side length of each blocking square.
  static const constexpr int kBlock = <?=$block?>;

  // The length of each column in the data matrix.
  static const constexpr int kHeight = <?=$size?>;

  // Whether the diagonal entries are kept.
  static const constexpr bool kDiag = <?=$diag?>;

 private:
  // The number of threads being used.
  int num_threads;

  // The mutex for fragments.
  mutex m_fragments;

<?  foreach ($states_ as $i => $state) { ?>
  const Mat<Type>& data_<?=$i?>;

  const Col<Key>& keys_<?=$i?>;

  Vec<Type> mean_<?=$i?>;

  Vec<Type> sdev_<?=$i?>
<?  } ?>

 public:
  <?=$className?>(<?=const_typed_ref_args($states)?>)
      : iteration(0),
<?  foreach (array_keys($states) as $i => $name) { ?>
        data_<?=$i?>(<?=$name?>.GetData()),
        keys_<?=$i?>(<?=$name?>.GetKeys()),
<?  } ?>
        m_fragments() {
  }

  void PrepareRound(WorkUnits& workers, int num_threads) {
  }

  void DoStep(Task& task, cGLA& gla) {
  }

  int GetNumFragments() {
<?  foreach ($states_ as $i => $state) { ?>
  mean_<?=$i?> = mean(data_<?=$i?>);
  sdev_<?=$i?> = stddev(data_<?=$i?>, 1);
<?  } ?>

<?  if ($single) { ?>
    int ratio = (data_1.n_cols - 1) / kBlock + 1;
    int n_fragments = ratio * (ratio + 1) / 2;
    return n_fragments;
<?  } else { ?>
    int ratio_1 = (data_1.n_cols - 1) / kBlock + 1;
    int ratio_2 = (data_2.n_cols - 1) / kBlock + 1;
    int n_fragments = ratio_1 * ratio_2;
    return n_fragments;
<?  } ?>
  }

  Iterator* Finalize(long fragment) {
<?  if ($single) { ?> {
    // The co-ordinates of the current block in the blocking grid are computed.
    int block_col = (sqrt(1 + 8 * fragment) - 1) / 2;
    int block_row = fragment - block_col * (block_col + 1) / 2;

    // The span of this block with regard to the covariance matrix.
    // Both intervas are closed on the left and open on the right.
    int first_col = kBlock * block_col;
    int first_row = kBlock * block_row;
    int final_col = std::min(count, first_col + kBlock);
    int final_row = std::min(count, first_row + kBlock);

    // The block, i.e. the covariance submatrix, is computed.
    Mat<Type> col_items = data.cols(first_col, final_col - 1);
    Mat<Type> row_items = data.cols(first_row, final_row - 1).t();

    Row<Type> col_means =   means.subvec(first_col, final_col - 1);
    Col<Type> row_means =   means.subvec(first_row, final_row - 1).t();
    Row<Type> col_stds  = stddevs.subvec(first_col, final_col - 1);
    Col<Type> row_stds  = stddevs.subvec(first_row, final_row - 1).t();

    Mat<Type> block = (row_items * col_items / kHeight - row_means * col_means)
                    / (row_stds * col_stds);

    // If the block lies on the main diagonal, not all of its elements are used.
    // Because the co-variance matrix is symmetric, only the upper triangular
    // portion of the matrix is returned.
    bool diagonal = (block_col == block_row);

    return new Iterator(first_col, first_row, n_cols, n_rows, diagonal,
                        std::move(block));
  }
<?  } else { ?>
    // The co-ordinates of the current block in the blocking grid are computed.
    int ratio = (data_1.n_cols - 1) / kBlock + 1;
    int block_col = fragment % ratio
    int block_row = fragment / ratio;

    // The span of this block with regard to the covariance matrix.
    // Both intervas are closed on the left and open on the right.
    int first_col = kBlock * block_col;
    int first_row = kBlock * block_row;
    int final_col = std::min(data_1.n_col, first_col + kBlock);
    int final_row = std::min(data_2.n_col, first_row + kBlock);

    // The block, i.e. the covariance submatrix, is computed.
    Mat<Type> col_items = data_1.cols(first_col, final_col - 1);
    Mat<Type> row_items = data_2.cols(first_row, final_row - 1).t();

    Row<Type> col_mean = mean_1.subvec(first_col, final_col - 1);
    Col<Type> row_mean = mean_2.subvec(first_row, final_row - 1).t();
    Row<Type> col_sdev = sdev_1.subvec(first_col, final_col - 1);
    Col<Type> row_sdev = sdev_2.subvec(first_row, final_row - 1).t();

    Mat<Type> block = (row_items * col_items / kHeight - row_means * col_means)
                    / (row_stds * col_stds);

    bool diagonal = false;

    return new Iterator(first_col, first_row, n_cols, n_rows, diagonal,
                        std::move(block));
<?  } ?>
  }

  bool GetNextResult(Iterator* it, <?=typed_ref_args($outputs_)?>) {
<?  if ($single) { ?>
    if (it->col == it->n_cols)
      return false;
    x = keys_1(it->row +it->row_shift);
    y = keys_1(it->col + it->col_shift);
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
<?  } else { ?>
    if (it->col == it->n_cols)
      return false;
    x = keys_1(it->row +it->row_shift);
    y = keys_2(it->col + it->col_shift);
    val = it->block(it->row, it->col);
    if (it->row == it->n_rows - 1) {
      it->row = 0;
      it->col++;
    } else {
      it->row++;
    }
    return true;
  }
<?  } ?>
};

<?
    return [
        'kind'            => 'GIST',
        'name'            => $className,
        'system_headers'  => $sys_headers,
        'user_headers'    => $user_headers,
        'lib_headers'     => $lib_headers,
        'libraries'       => $libraries,
        'extra'           => $extra,
        'iterable'        => true,
        'intermediate'    => false,
        'output'          => $outputs,
        'result_type'     => $result_type,
    ];
}
?>
