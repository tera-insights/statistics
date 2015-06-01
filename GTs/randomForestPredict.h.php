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
    // Class name randomly generated.
    $className = generate_name("RFP");

    // Initialization of local variables from template arguments.
    $file = get_default($t_args, 'file', false);

    grokit_assert($file || count($states) > 0,
                 "The model must be passed in either the state or a file.");

    // Naming the inputs.
    $inputs_ = array_combine(['x'], $inputs);

    $type = array_get_index($states, 0)->get('type');
    array_set_index($outputs, 0, $type);
    $outputs_ = array_combine(['y'], $outputs);

    $sys_headers  = ['armadillo', 'opencv/cv.h', 'opencv/ml.h'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo', 'opencv_core', 'opencv_ml'];
?>

using namespace cv;
using namespace std;

class <?=$className?>;

<?  $constantState = lookupResource(
        'statistics::Random_Forest_Predict_Constant_State',
        ['className' => $className, 'states' => $states, 'file' => $file]
    ); ?>

class <?=$className?> {
 private:
  // The constant state containing the forest.
  const <?=$constantState?>& constant_state;

  // The vector to be filled with the predictors.
  Mat sample;

 public:
  <?=$className?>(const <?=$constantState?>& state)
      : constant_state(state) {
  }

  // The sample is most likely being placed directly on top of the input's
  // memory, which could be dangerious.
  bool ProcessTuple(<?=process_tuple_args($inputs_, $outputs_)?>) {
    sample = Mat(x.n_cols, x.n_rows, CV_32F, (void*) x.memptr());
    y = constant_state.forest.predict(sample, Mat());
    return true;
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
        'iterable'        => false,
        'input'           => $inputs,
        'output'          => $outputs,
        'result_type'     => 'single',
    ];
}
?>
