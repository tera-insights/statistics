<?
function Random_Forest(array $t_args, array $inputs, array $outputs)
{
    // Class name is randomly generated.
    $className = generate_name("RF");

    // Initializiation of argument names.
    $inputs_ = array_combine(['x', 'y'], $inputs);

    // Initialization of local variables from template arguments
    $scale = get_default($t_args, 'scale', 2);
    $width = get_default($t_args, 'length', 100);

    $file          = get_default($t_args, 'file',           false);
    $maxDepth      = get_default($t_args, 'max.depth',      25);
    $sampleCount   = get_default($t_args, 'min.sample',     100);
    $nodeEpsilon   = get_default($t_args, 'node.epsilon',   0.01);
    $maxCategories = get_default($t_args, 'max.categories', 15);
    $numVars       = get_default($t_args, 'num.vars',       0);
    $numTrees      = get_default($t_args, 'num.trees',      0);
    $treeEpsilon   = get_default($t_args, 'tree.epsilon',   0);

    grokit_assert($numTrees + $treeEpsilon > 0,
                  'Random Forest: no stopping criterion given.');

    if ($numTrees > 0)
        $stopping = 'CV_TERMCRIT_ITER';
        if ($treeEpsilon > 0)
            $stopping .= ' | CV_TERMCRIT_EPS';
    else
        $stopping = 'CV_TERMCRIT_EPS';

    $vector = $inputs_['x'];
    $height = $vector->get('size');
    $types = array_values(array_merge($vector->get('inputs'), [$inputs_['y']]));
    foreach ($types as &$type) {
        grokit_assert($type->is('categorical') || $type->is('numeric'),
                      'Random Forest: unable to use type $type');
        $type = $type->is('numeric') ? 'CV_VAR_NUMERICAL' : 'CV_VAR_CATEGORICAL';
    }

    $sys_headers  = ['armadillo', 'opencv/cv.h', 'opencv/ml.h'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo', 'opencv_core', 'opencv_ml'];
    $extra        = ['type' => $inputs_['y']];
    $result_type  = ['state'];
?>

using namespace arma;
using namespace std;
using Mat = cv::Mat;

class <?=$className?>;

class <?=$className?> {
 public:

  // The length of each column in the data matrix.
  static const constexpr unsigned int kHeight = <?=$height?>;

  // The initial width of the data matrix.
  static const constexpr unsigned int kWidth = <?=$width?>;

  // The proportion at which the dynamic matrix grows.
  static const constexpr unsigned int kScale = <?=$scale?>;

 private:
  // The data matrix being constructed item by item. The width of this matrix
  // is increased when necessary as per a dynamic array.
  fmat features;

  // The vector of responses that correspond to each column in features.
  fvec response;

  // The number of rows processed by this state.
  unsigned int count;

  // The model to be trained.
  CvRTrees forest;

 public:
  <?=$className?>()
      : features(kHeight, kWidth),
        response(kWidth),
        count(0),
        forest() {
  }

  // Basic dynamic array allocation.
  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    if (count == features.n_cols) {
      features.resize(kHeight, kScale * features.n_cols);
      response.resize(kScale * response.n_elem);
    }
    features.col(count) = fvec(x);
    response(count) = y;
    count++;
  }

  // Empty rows are stripped such that white space will only ever be at the end
  // of both response and features.
  void AddState(<?=$className?> &other) {
    features.resize(kHeight, count);
    response.resize(count);
    features = join_rows(features, other.features);
    response = join_rows(response, other.response);
    count += other.count;
  }

  void FinalizeState() {
    // The remaining whitespace is stripped.
    features.resize(kHeight, count);
    response.resize(count);
    wall_clock timer;
    timer.tic();
    cout << "beginning training" << endl;
    cout << "features: " << features.n_rows << " by " << features.n_cols << endl;
    cout << "response: " << response.n_rows << " by " << response.n_cols << endl;
    Mat trainData = Mat(features.n_cols, features.n_rows, CV_32F, features.memptr());
    Mat responses = Mat(response.n_cols, response.n_rows, CV_32F, response.memptr());
    CvRTParams params = CvRTParams(
        <?=$maxDepth?>, <?=$sampleCount?>, <?=$nodeEpsilon?>, false,
        <?=$maxCategories?>, 0, false, <?=$numVars?>, <?=$numTrees?>,
        <?=$treeEpsilon?>, <?=$stopping?>
    );
    Mat types = Mat(kHeight + 1, 1, CV_8U);
<?  foreach ($types as $counter => $type) { ?>
    types.at<uchar>(<?=$counter?>, 0) = <?=$type?>;
<?  } ?>
    forest.train(trainData, CV_ROW_SAMPLE, responses,
                 Mat(), Mat(), types, Mat(), params);
    cout << "finished training in " << timer.toc() << " seconds." << endl;
    cout << "var_type: " << Mat(forest.get_tree(0)->get_data()->var_type) << endl;
    cout << "var_idx: " << Mat(forest.get_tree(0)->get_data()->var_idx) << endl;
    cout << "cat_map: " << Mat(forest.get_tree(0)->get_data()->cat_map) << endl;
    cout << "cat_ofs: " << Mat(forest.get_tree(0)->get_data()->cat_ofs) << endl;
    cout << "cat_count: " << Mat(forest.get_tree(0)->get_data()->cat_count) << endl;
<?  if ($file) { ?>
    FileStorage file = FileStorage("<?=$file?>", FileStorage::WRITE);
    forest.write(*file, "rtree");
    file.release();
<?  } ?>
  }

 public:
  const CvRTrees& GetForest() const {
    return forest;
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
        'extra'             => $extra,
        'iterable'          => false,
        'input'             => $inputs,
        'output'            => $outputs,
        'finalize_as_state' => true,
        'result_type'       => $result_type,
    ];
}
?>
