<?
function Random_Forest_Predict_Constant_State($t_args)
{
    // Initialization of local variables from template arguments.
    $className = $t_args['className'];
    $states    = $t_args['states'];
    $file      = $t_args['file'];

    if (!$file)
        $states_ = array_combine(['state'], $states);

    // Return values.
    $sys_headers  = ['armadillo', 'opencv/cv.h', 'opencv/ml.h'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo', 'opencv_core', 'opencv_ml'];
?>

using namespace cv;

class <?=$className?>ConstantState {
 public:
  friend class <?=$className?>;

 private:
  // The forest used in prediction.
  CvRTrees forest;

 public:
<?  if ($file) { ?>
  <?=$className?>ConstantState() {
    FileStorage file = FileStorage("<?=$file?>", FileStorage::READ);
    forest.read(*file, *file.getFirstTopLevelNode());
  }
<?  } else { ?>
  <?=$className?>ConstantState(<?=const_typed_ref_args($states_)?>)
      : forest(state.GetForest()) {
  }
<?  } ?>
};
<?
    return [
        'kind'           => 'RESOURCE',
        'name'           => $className . 'ConstantState',
        'system_headers' => $sys_headers,
        'user_headers'   => $user_headers,
        'lib_headers'    => $lib_headers,
        'libraries'      => $libraries,
        'configurable'   => false,
    ];
}

//  Copyright 2014 Tera Insights, LLC. All Rights Reserved.
function Random_Forest_Predict($t_args, $inputs, $outputs, $states)
{
    // Class name is randomly generated.
    $className = generate_name("RFP");

    // Initialization of local variables from template arguments.
    $file  = get_default($t_args, 'file',   false);
    $scale = get_default($t_args, 'scale',  2);
    $width = get_default($t_args, 'length', 100);

    grokit_assert($file || count($states) > 0,
                 "The model must be passed in either the state or a file.");

    // Naming the inputs.
    $inputs_ = array_combine(['x'], $inputs);

    $output = array_get_index($states, 0)->get('type');
    $regression = $output->is('numeric');
    if (!$regression)
      $cardinality = $output->get('cardinality');
    else
      $cardinality = 1;
    array_set_index($outputs, 0, $output);
    $outputs_ = array_combine(['y'], $outputs);

    $sys_headers  = ['armadillo', 'opencv/cv.h', 'opencv/ml.h'];
    $user_headers = [];
    $lib_headers  = ['tree.h'];
    $libraries    = ['armadillo', 'opencv_core', 'opencv_ml'];
?>

using namespace std;

class <?=$className?>;

<?  $constantState = lookupResource(
        'statistics::Random_Forest_Predict_Constant_State',
        ['className' => $className, 'states' => $states, 'file' => $file]
    ); ?>

class <?=$className?> {
 public:
  // The type of each tree.
  using Tree = GiDTree;

  // The type of matrices being trained and predicted on.
  using Mat = cv::Mat;

  // The initial width of the results matrix.
  static const constexpr int kWidth = <?=$width?>;

  // The proportion at which the dynamic matrix grows.
  static const constexpr int kScale = <?=$scale?>;

<?  if (!$regression) { ?>
  // The cardinality of the output.
  static const constexpr int kCardinality = <?=$cardinality?>;
<?  } ?>

 private:
  // The constant state containing the forest.
  const <?=$constantState?>& constant_state;

  // The number of trees in the random forest.
  const int num_trees;

<?  if ($regression) { ?>
  // The sum of the predictions for each user.
<?  } else { ?>
  // The distribution of votes per category for each user.
<?  } ?>
  mat total;

  // The current iteration of the GT for the current chunk.
  int iteration;

  // The index of the current item for this chunk.
  int index;

  // The vector to be filled with the predictors.
  Mat sample;

  // The tree currently be used for prediction.
  Tree tree;

 public:
  <?=$className?>(const <?=$constantState?>& state)
      : constant_state(state),
        num_trees(constant_state.forest.get_tree_count()),
        total(kCardinality, count) {
  }

  void StartChunk() {
    iteration = 0;
    ResetChunk();
  }

  // The sample is most likely being placed directly on top of the input's
  // memory, which could be dangerious.
  bool ProcessTuple(<?=process_tuple_args($inputs_, $outputs_)?>) {
    // First we increase the size of the results if necessary.
    if (index > total.n_cols)
        total.resize(total.n_rows, total.n_cols * kScale);
    // The prediction is performed and stored in the result.
    sample = Mat(x.n_cols, x.n_rows, CV_32F, (void*) x.memptr());
    y = tree.predict(x);
    // The respective total is increased.
<?  if ($regression) { ?>
    total(0, index) += y;
<?  } else { ?>
    total(y, index)++;
<?  } ?>
    // On the final iteration, the result is returned.
    if (iteration == num_trees) {
<?  if ($regression) { ?>
      y = total(0, index) / num_trees;
<?  } else { ?>
      total.col(index).max(y);
<?  } ?>
      return true;
    }
    // On all other iterations, no result is returned.
    return false;
  }

  bool ShouldIterate() {
    iteration++;
    if (iteration == num_trees)
      return false;
    ResetChunk();
    return true;
  }

 private:
  // This is used whenever a new iteration is begun.
  ResetChunk() {
    index = 0;
    tree = Tree(constant_state.forest.get_tree(iteration));
  }
};

<?
    return [
        'kind'            => 'GT',
        'name'            => $className,
        'generated_state' => $constantState,
        'system_headers'  => $sys_headers,
        'user_headers'    => $user_headers,
        'lib_headers'     => $lib_headers,
        'libraries'       => $libraries,
        'iterable'        => true,
        'input'           => $inputs,
        'output'          => $outputs,
        'result_type'     => 'single',
    ];
}
?>
