<?
function Random_Forest_Batch($t_args, $outputs, $states)
{
    // Class name randomly generated.
    $className = generate_name('RFB');

    $states_ = array_combine(['training', 'predicting'], $states);
    $output = array_values($states_['predicting']->input());
    $input  = array_values($states_['training']->input());

    $tuple = array_get_index($states_['predicting']->input(), 1);

    // Initialization of local variables from template arguments
    $file          = get_default($t_args, 'file',           false);
    $maxDepth      = get_default($t_args, 'max.depth',      25);
    $sampleCount   = get_default($t_args, 'min.sample',     100);
    $nodeEpsilon   = get_default($t_args, 'node.epsilon',   0.01);
    $maxCategories = get_default($t_args, 'max.categories', 15);
    $numVars       = get_default($t_args, 'num.vars',       0);

    $numTrees = $t_args['num.trees'];

    // Setting output types.
    $outputs_ = [];

    $types = array_values($output[1]->get('types'));
    $length = count($types);
    for ($i = 0; $i < $length; $i++)
        $outputs_["extra$i"] = $types[$i];

    $types = array_values($output[0]->get('inputs'));
    $height = count($types);
    for ($i = 0; $i < $height; $i++)
        $outputs_["input$i"] = $types[$i];

    $outputs_['output'] = $output = array_get_index($input[1]->get('types'), 0);
    $outputs = array_combine(array_keys($outputs), $outputs_);

    $types = array_values(array_merge($input[0]->get('inputs'), [$output]));
    foreach ($types as &$type) {
        grokit_assert($type->is('categorical') || $type->is('numeric'),
                      'Random Forest: unable to use type $type');
        $type = $type->is('numeric') ? 'CV_VAR_NUMERICAL' : 'CV_VAR_CATEGORICAL';
    }

    $regression = $output->is('numeric');
    if (!$regression)
      $cardinality = $output->get('cardinality');

    $sys_headers  = ['armadillo', 'vector', 'mutex'];
    $user_headers = [];
    $lib_headers  = ['tree.h'];
    $libraries    = ['armadillo', 'opencv_core', 'opencv_ml'];
    $extra        = [];
    $result_type  = ['fragment'];
?>

using namespace cv;
using namespace arma;
using namespace std;

class AnswerGLA {
 private:
  // Number of iterations finished.
  bool answer;

 public:
  AnswerGLA(bool answer)
      : answer(answer) {
  }

  void AddState(AnswerGLA other) {};

  bool ShouldIterate() {
    return answer;
  }
};

class <?=$className?>;

class <?=$className?> {
 public:
  struct Iterator {
    // The fragment ID for this tree, corresponding to the tree index.
    long fragment;

    // The index for the current output of this fragment.
    long index;
  };

  struct Task {
    // The index corresponding to that of the scheduler allocating this task.
    int index;
  };

  struct LocalScheduler {
    // The thread index of this scheduler.
    int index;

    // Whether this scheduler has scheduled its single task.
    bool finished;

    LocalScheduler(int index)
        : index(index),
          finished(false) {
    }

    bool GetNextTask(Task& task) {
      bool ret = !finished;
      printf("Getting task from scheduler %d: %d\n", index, ret);
      task.index = index;
      finished = true;
      return ret;
    }
  };

  // The type of each tree.
  using Tree = GiDTree;

  // The type of matrices being trained and predicted on.
  using Mat = cv::Mat;

  // The type of the extra tuples.
  using Tuple = <?=$tuple?>;

  // The inner GLA being used.
  using cGLA = AnswerGLA;

  // The type of the workers.
  using WorkUnit = pair<LocalScheduler*, cGLA*>;

  // The type of the container for the workers.
  using WorkUnits = vector<WorkUnit>;

  // The length of each column in the items matrix.
  static const constexpr int kHeight = <?=$height?>;

  // The length of each tuple in the extra vector.
  static const constexpr int kLength = <?=$length?>;

  // The number of trees in the forest.
  static const constexpr int kNumTrees = <?=$numTrees?>;

<?  if (!$regression) { ?>
  // The cardinality of the output.
  static const constexpr int kCardinality = <?=$cardinality?>;
<?  } ?>

 private:
  // The current iteration of the GIST.
  int iteration;

  // The number of items being processed;
  long count;

  // The matrix containing the data used for training.
  const fmat& training;

  // The vector of responses corresponding to the training data.
  fvec response;

  // The matrix containing the data used for predicting.
  const mat& predicting;

  // The vector containing the extra attributes to pass through when predicting.
  const vector<Tuple>& extra;

  // The random forest used to do prediction, created by combining the forests
  // created in each work unit.
  array<CvDTree*, kNumTrees> forest;

  // The forests trained locally. This is used to delete them at the end.
  vector<CvRTrees*> forests;

  // The number of threads being used.
  int num_threads;

<?  if ($regression) { ?>
  // The sum of the predictions for each user.
  vec total;
<?  } else { ?>
  // The distribution of votes per category for each user.
  mat total;
<?  } ?>

  // The mutex for fragments.
  mutex m_fragments;

 public:
  <?=$className?>(<?=const_typed_ref_args($states_)?>)
      : iteration(0),
        count(predicting.GetCount()),
        training(training.GetMatrix()),
        response((float*) training.GetTuples().data(), training.GetCount()),
        predicting(predicting.GetMatrix()),
        extra(predicting.GetTuples()),
        total(<?=$regression ? '' : 'kCardinality, '?>count),
        m_fragments() {
    cout << "constructed gist state" << endl;
    cout << "count " << predicting.GetCount() << endl;
    cout << "training " << this->training.n_rows << " x " << this->training.n_cols << endl;
    cout << "response " << this->response.n_rows << " x " << this->response.n_cols << endl;
    cout << "predicting " << this->predicting.n_rows << " x " << this->predicting.n_cols << endl;
  }

  ~<?=$className?>() {
    for (auto forest : forests)
      delete forest;
  }

  void PrepareRound(WorkUnits& workers, int num_threads) {
    this->num_threads = num_threads;
    int num_workers = iteration ? kNumTrees : num_threads;
    // if (iteration == 0)
    //   num_workers = this->num_threads = 2;
    cout << "Beginning round " << iteration << " with " << num_workers << " workers." << endl;
    for (int counter = 0; counter < num_workers; counter++)
      workers.push_back(WorkUnit(new LocalScheduler(counter), new cGLA(!iteration)));
    iteration++;
  }

  void DoStep(Task& task, cGLA& gla) {
    // printf("doing step during iteration %d\n", iteration);
    if (iteration == 1) {
      int num_trees = (task.index + 1) * kNumTrees / num_threads
                    - task.index * kNumTrees / num_threads;
      // printf("worker %d is training %d num trees\n", task.index, num_trees);
      CvRTParams params(<?=$maxDepth?>, <?=$sampleCount?>, <?=$nodeEpsilon?>,
                        false, <?=$maxCategories?>, 0, false, <?=$numVars?>,
                        num_trees, 0, CV_TERMCRIT_ITER);
      Mat trainData = Mat(training.n_cols, training.n_rows, CV_32F, (void*) training.memptr());
      Mat responses = Mat(response.n_cols, response.n_rows, CV_32F, response.memptr());
      Mat types = Mat(kHeight + 1, 1, CV_8U);
<?  foreach ($types as $counter => $type) { ?>
      types.at<uchar>(<?=$counter?>, 0) = <?=$type?>;
<?  } ?>
      CvRTrees* forest = new CvRTrees();
      forest->train(trainData, CV_ROW_SAMPLE, responses,
                    Mat(), Mat(), types, Mat(), params);
      int index = task.index * kNumTrees / num_threads;
      printf("task %d trainined %d trees: %d to %d\n",
             task.index, num_trees, index, index + num_trees - 1);
      for (int counter = 0; counter < num_trees; counter++)
        this->forest[index + counter] = forest->get_tree(counter);
      { // Locking on field variables
        unique_lock<mutex> guard(m_fragments);
        forests.push_back(forest);
      }
    } else {
      Tree tree(forest[task.index]);
<?  if ($regression) { ?>
      vec predictions(count);
      for (int index = 0; index < count; index++)
        predictions(index) = tree.predict(predicting.col(index));
<?  } else { ?>
      mat predictions(cardinality, count);
      for (int index = 0; index < count; index++)
        predictions(tree.predict(predicting.col(index)), index)++;
<?  } ?>
      { // Locking on field variables
        unique_lock<mutex> guard(m_fragments);
        total += predictions;
      }
    }
  }

  int GetNumFragments() {
<?  if ($regression) { ?>
    total /= kNumTrees;
<?  } ?>
    return 25 * (iteration - 1);
  }

  Iterator* Finalize(long fragment) {
    return new Iterator{fragment, fragment * count / 25};
  }

  bool GetNextResult(Iterator* it, <?=typed_ref_args($outputs_)?>) {
    if (it->index == (it->fragment + 1) * count / 25)
      return false;
<?  for ($i = 0; $i < $length; $i++) { ?>
    extra<?=$i?> = get<<?=$i?>>(extra[it->index]);
<?  } ?>
<?  for ($i = 0; $i < $height; $i++) { ?>
    input<?=$i?> = predicting(<?=$i?>, it->index);
<?  } ?>
<?  if ($regression) { ?>
    output = total(it->index);
<?  } else { ?>
    total.col(it->index).max(output);
<?  } ?>
    it->index++;
    return true;
  }
};

typedef <?=$className?>::Iterator <?=$className?>_Iterator;

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
