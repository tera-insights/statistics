<?
function Random_Forest_Batch_Constant_State($t_args)
{
    // Initialization of local variables from template arguments.
    $className = $t_args['className'];
    $states    = $t_args['states'];

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
  const CvRTrees& forest;

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

function Random_Forest_Batch($t_args, $inputs, $outputs, $states)
{
    // Class name randomly generated.
    $className = generate_name("RFB");

    // Whether there are extra attributes to pass through.
    $hasExtra = count($inputs) > 1;

    // Initialization of local variables from input names.
    if ($hasExtra)
        $inputs_ = array_combine(['extras', 'vector'], $inputs);
    else
        $inputs_ = array_combine(['vector'], $inputs);

    // Initialization of local variables from template arguments.
    $scale = get_default($t_args, 'scale',  2);
    $width = get_default($t_args, 'length', 100);
    $file  = get_default($t_args, 'file',   false);

    if ($hasExtra)
        $length = $inputs_['extras']->get('size');
    else
        $length = 0;
    $height = $inputs_['vector']->get('size');

    // Setting output types.
    $outputs_ = [];
    if ($hasExtra) {
        $types = array_values($inputs_['extras']->get('inputs'));
        for ($i = 0; $i < $length; $i++)
            $outputs_["extra$i"] = $types[$i];
    }
    $types = array_values($inputs_['vector']->get('inputs'));
    for ($i = 0; $i < $height; $i++)
        $outputs_["input$i"] = $types[$i];
    $outputs_["output"] = array_get_index($states, 0)->get('type');
    $outputs = array_combine(array_keys($outputs), $outputs_);

    $regression = $outputs_['output']->is('numeric');
    if (!$regression)
      $cardinality = $outputs_['output']->get('cardinality');

    $sys_headers  = ['armadillo', 'limits', 'mutex'];
    $user_headers = [];
    $lib_headers  = ['tree.h'];
    $libraries    = ['armadillo'];
    $extra        = [];
    $result_type  = ['fragment'];
?>

using namespace arma;
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

  // The length of each column in the items matrix.
  static const constexpr unsigned int kHeight = <?=$height?>;

  // The length of each column in the extra matrix.
  static const constexpr unsigned int kLength = <?=$length?>;

  // The initial width of the data matrix.
  static const constexpr unsigned int kWidth = <?=$width?>;

  // The proportion at which the dynamic matrix grows.
  static const constexpr unsigned int kScale = <?=$scale?>;

<?  if (!$regression) { ?>
  // The cardinality of the output.
  static const constexpr unsigned int kCardinality = <?=$cardinality?>;
<?  } ?>

  struct Iterator {
    // The fragment ID for this tree, corresponding to the tree index.
    unsigned int fragment;

    // Whether this fragment has already done its predictions.
    bool predicted;

    Iterator(int fragment)
        : fragment(fragment),
          predicted(false) {
    }
  };

 private:
  // The constant state containing the forest.
  const <?=$constantState?>& constant_state;

  // The matrices containing the predictors and passed through attributes.
  mat items, extra;

  // The number of inputs processed by this state.
  unsigned int count;

  // The number of fragments remaining.
  unsigned int fragments_remaining;

<?  if ($regression) ?>
  // The sum of the predictions for each user.
  vec total;
<?  } else { ?>
  // The distribution of votes per category for each user.
  mat total;
<?  } ?>

  // The mutex for fragments.
  mutex m_fragments;

  // The index of the current output.
  unsigned int index;

 public:
  <?=$className?>(const <?=$constantState?>& state)
      : constant_state(state),
        items(kHeight, kWidth),
        extra(kLength, kWidth),
        count(0),
        m_fragments(),
        index(0) {
  }

  // Basic dynamic array allocation.
  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    if (count == items.n_cols) {
      items.resize(kHeight, kScale * items.n_cols);
<?  if ($hasExtra) { ?>
      extra.resize(kLength, kScale * extra.n_elem);
<?  } ?>
    }
    items.col(count) = vector;
<?  if ($hasExtra) { ?>
    extra.col(count) = extras;
<?  } ?>
    count++;
  }

  // Empty rows are stripped such that white space will only ever be at the end
  // of both items and extra;
  void AddState(<?=$className?> &other) {
    items.resize(kHeight, count);
    items.insert_cols(count, other.items);
<?  if ($hasExtra) { ?>
    extra.resize(kLength, count);
    extra.insert_cols(count, other.extra);
<?  } ?>
    count += other.count;
  }

  // There is one fragment per tree.
  int GetNumFragments() {
    // The remaining whitespace is stripped.
    items.resize(kHeight, count);
<?  if ($hasExtra) { ?>
    extra.resize(kLength, count);
<?  } ?>

<?  if ($regression) { ?>
    total.set_size(count);
<?  } else { ?>
    votes.set_size(cardinality, count);
<?  } ?>
    total.fill(0);

    return fragments_remaining = constant_state.forest.get_tree_count();
  }

  Iterator* Finalize(int fragment) {
    return new Iterator(fragment);
  }

  bool GetNextResult(Iterator* it, <?=typed_ref_args($outputs_)?>) {
    if (!it->predicted) {
      Tree tree (constant_state.forest.get_tree(it->fragment));
<?  if ($regression) { ?>
      vec predictions (count);
      for (int index = 0; index < count; index++)
        predictions(index) = tree.predict(items.col(index));
<?  } else { ?>
      mat predictions (cardinality, count);
      for (int index = 0; index < count; index++)
        predictions(tree.predict(items.col(index)), index)++;
<?  } ?>

      { // Locking on field variables
        unique_lock<mutex> guard (m_fragments);
        total += predictions;
        fragments_remaining--;

        if (fragments_remaining)
          return false;
      }

<?  if ($regression) { ?>
      total /= constant_state.forest.get_tree_count();
<?  } ?>
      it->predicted = true;
    }

    if (index == count)
      return false;

<?  for ($i = 0; $i < $length; $i++) { ?>
    extra<?=$i?> = extra(<?=$i?>, index);
<?  } ?>
<?  for ($i = 0; $i < $height; $i++) { ?>
    input<?=$i?> = items(<?=$i?>, index);
<?  } ?>
<?  if ($regression) ?>
    output = total(index);
<?  } else { ?>
    total.col(index).max(output);
<?  } ?>

    index++;
    return true;
  }
};

typedef <?=$className?>::Iterator <?=$className?>_Iterator;

<?
    return [
        'kind'            => 'GLA',
        'name'            => $className,
        'generated_state' => $constantState,
        'system_headers'  => $sys_headers,
        'user_headers'    => $user_headers,
        'lib_headers'     => $lib_headers,
        'libraries'       => $libraries,
        'extra'           => $extra,
        'iterable'        => false,
        'input'           => $inputs,
        'output'          => $outputs,
        'result_type'     => $result_type,
    ];
}
?>
