<?
function G_Means_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className  = $t_args['className'];
    $dimension  = $t_args['dimension'];
?>
#define ROOT_HALF_PI (boost::math::double_constants::root_half_pi)

using namespace arma;

// These states are detailed in the cluster class.
enum class <?=$className?>State {
  kMeans,
  gMeans,
  kMeansChildren,
  jarqueBera,
  activeParent,
  noChildren,
  inactiveParent,
};

class <?=$className?>;

// This is a constant node. It does not maintain any statistics that will change
// during an iteration. It is instead copied over into the Cluster class.

// The variables here lack description. They are analagous to those in the cluster.

class <?=$className?>Node {
 typedef <?=$className?>State State;

 private:
  int count;
  vec center;
  <?=$className?>State state;
  <?=$className?>Node *child1;
  <?=$className?>Node *child2;

 public:
  friend class <?=$className?>Cluster;
  friend class <?=$className?>;

  <?=$className?>Node(const vec &_center)
      : count(0),
        center(_center),
        state(State::kMeans),
        child1(nullptr),
        child2(nullptr) {
  }

  // A typical deconstructor for the allocated children.
  ~<?=$className?>Node() {
    if (child1 != nullptr) {
      delete child1;
      child1 = nullptr;
    }
    if (child2 != nullptr) {
      delete child2;
      child2 = nullptr;
    }
  }

  // This is called at the end of every iteration and decides if a parent has
  // become inactive.
  bool MakeInactive() {
    if (state == State::noChildren || state == State::inactiveParent) {
      return true;
    } else if (state == State::activeParent && child1->MakeInactive() && child2->MakeInactive()) {
      state = State::inactiveParent;
      return true;
    } else {
      return false;
    };
  }

  // TODO: Change this to JSON
  Json::Value GetOutput(int level) const {
    Json::Value data(Json::objectValue);
    //for (int counter = 0; counter < level; counter ++)
      //cout << "    ";
    switch (state) {
      case State::kMeans:
        data["category"] = "kmeans";
        break;
      case State::gMeans:
        data["category"] = "gMeans";
        break;
      case State::kMeansChildren:
        data["category"] = "kMeansChildren";
        break;
      case State::jarqueBera:
        data["category"] = "jarqueBera";
        break;
      case State::activeParent:
        data["category"] = "activeParent";
        break;
      case State::inactiveParent:
        data["category"] = "inactiveParent";
        break;
      case State::noChildren:
        data["category"] = "noChildren";
        break;
    }
    //cout << center << endl;
<?  foreach (range(0, $dimension - 1) as $counter) { ?>
    data["center"].append(center(<?=$counter?>));
<?  }  ?>
    if (child1 != nullptr) data["offspring1"] = child1->GetOutput(level + 1);
    if (child2 != nullptr) data["offspring2"] = child2->GetOutput(level + 1);
    data["depth"] = level;
    data["count"] = count;
    return data;
  }
};

class <?=$className?>ConstantState{
 private:
  // TODO: Maybe add a joint iteration. Probably too complicated for the average
  // user because the intermixing of splitting, kmeans, and computations.
  // The current iteration which is to be compared to $axIteration
  long iteration;

  // The current depth of the tree which is to be compared to $maxDepth.
  long depth;

  // The root of the tree representing the hierarchy of clusters.
  <?=$className?>Node root;

 public:
  friend class <?=$className?>;

  // The root is initialized with a zero vector and shifts to the mean during
  // the first iteration.
  <?=$className?>ConstantState()
      : iteration(0),
        depth(0),
        root(zeros<vec>(<?=$dimension?>)) {
  }
};
<?
  return array(
        'kind'           => 'RESOURCE',
        'name'           => $className . 'ConstantState',
        'system_headers' => array('armadillo', 'vector'),
        'user_headers'   => array('json.h'),
        'user_headers'   => array(),
    );
}

function G_Means(array $t_args, array $inputs, array $outputs)
{
    // Setting output type
    array_set_index($outputs, 0, lookupType("base::JSON"));

    // Class name is randomly generated
    $className = generate_name("HGM");

    // The dimension of the data, i.e. how many elements are in each item.
    $dimension = count($inputs);

    // Processing of template arguments

    // The parameter p for the Minkowski distance metric used. It is not
    // technically a metric is p < 1 but this is allowable.
    $pMinkowski = $t_args['p'];

    // The maximum allowed change allowed for convergence; It should be 0 for
    // convergence only when the clusters are stationary.
    $epsilon = $t_args['epsilon'];

    // The significance level of the Jarque-Bera test.
    $alpha = $t_args['alpha'];

    $codeArray = array_combine(array_keys($inputs), range(0, $dimension - 1));
?>
using namespace arma;
using namespace std;

class <?=$className?>Cluster {
 typedef <?=$className?>State State;

 private:
  // This records the current state of the algorithm to be used in a switch.
  // 0: Run k-means. Multiple iterations until convergence. Calculates v.
  // 1: Calculate XTX and split for G-Means. Single iteration.
  // 2: Waiting for children k-means to converge.
  // 3: Transform the data with v and compute the mean.
  // 4: Calculate the Jarque Bera statistic.
  // 5: The above is complete and this cluster has children.
  // 6: This cluster has been completed and has no children.
  // 7: This cluster has children but all of its offspring are complete.
  // TODO: Maybe switch this to a set of booleans to gain branch prediction
  // instead of lookup table. Gets predicted almost perfectly due to low change.
  // TODO: Maybe merge 3 and 4 - empirically test this, especially because less projections.
  // It should be noted that this state is a bit weird because this cluster can
  // either be a parent or a child. The child is initialized by the parent with
  // state 0 and k-means is run on it. The child is then untouched until the
  // conclusion of G-Means on its parent, after which it is either deleted if
  // the parent is state 6 or kept and set to state 1 if the parent is state 5.
  // The child is then at state 1 and begins its own G-Means. The state is
  // switched to 7 during a check in ShouldIterate whichs checks if all of its
  // offspring are state 6 or 7.
  State state;

  // The center of this cluster.
  vec center;

  // The children nodes that partition this cluster.
  <?=$className?>Cluster *child1;
  <?=$className?>Cluster *child2;

  // The sum of the element in this cluster; used to calculate the centroid.
  vec sum;

  // The total number of elements in the cluster. Incremented by the AddItem
  // class function and used to calculate the centroid.
  long count;

  // This represents X^T * X where X is the data matrix for this cluster. It is
  // needed for principal component analysis in order to split the cluster.
  // It is the first thing done when a new cluster is made.
  mat XTX;

  // This represents a major axis of the cluster and is the difference between
  // the two child cluster centers. It has the same name as in the G-means paper
  vec v;

  // This is the norm of v, stored to avoid unnecessary computation.
  double v_norm;

  // These are used in the calculation of the Jarque-Bera statistics, i.e. the
  // skewness and kurtosis.

  // The following 4 variables stored the sum of the projections raised to
  // various values, as specified by their suffix.
  long double sum_first;
  long double sum_second;
  long double sum_third;
  long double sum_fourth;

  // The scalar value of the projection of the item unto the axis, v. It is
  // stored to avoid recomputation for each of the 3 above statistics.
  double projection;

  // This is used to aid calculation of the moments by storing intermediate
  // products.
  double dummy_term;

 public:
  friend class <?=$className?>;

  <?=$className?>Cluster(const <?=$className?>Node &node)
      : sum_first(0),
        sum_second(0),
        sum_third(0),
        sum_fourth(0),
        count(node.count),
        state(node.state),
        center(node.center),
        sum(zeros<vec>(<?=$dimension?>)),
        XTX (zeros<mat>(<?=$dimension?>, <?=$dimension?>)) {
    switch (state) {
      case State::kMeansChildren:
        child1 = new <?=$className?>Cluster(*node.child1);
        child2 = new <?=$className?>Cluster(*node.child2);
        break;
      case State::jarqueBera:
        child1 = new <?=$className?>Cluster(*node.child1);
        child2 = new <?=$className?>Cluster(*node.child2);
        v = child1->center - child2->center;
        v_norm = norm(v, <?=$pMinkowski?>);
        break;
      case State::activeParent:
      case State::inactiveParent:
        child1 = new <?=$className?>Cluster(*node.child1);
        child2 = new <?=$className?>Cluster(*node.child2);
        break;
      default:
        child1 = nullptr;
        child2 = nullptr;
    }
  }

  ~<?=$className?>Cluster() {
    if (child1 != nullptr) {
      delete child1;
      child1 = nullptr;
    }
    if (child2 != nullptr) {
      delete child2;
      child2 = nullptr;
    }
  }

  <?=$className?>Cluster(const <?=$className?>Cluster &)=delete;
  <?=$className?>Cluster & operator =(const <?=$className?>Cluster &)=delete;

  // This function recalculates the statistics of the cluster upon adding a
  // point, i.e. the data item. It is call AddItem only for sake of consistency
  // and this name is not relevant to Datapath - this function is only ever
  // called by AddItem of the GLA below.
  // TODO: Maybe change from switch to if-else block for prediction?
  void AddItem(vec &point) {
    switch (state) {
      case State::kMeans:
        count ++;
        sum += point;
        break;
      case State::gMeans:  // Add to the PVA matrix.
        XTX += (point - center) * trans(point - center);
        break;
      case State::kMeansChildren:  // Filter the data for k-means.
        if (child1->GetDistance(point) < child2->GetDistance(point))
          child1->AddItem(point);
        else
          child2->AddItem(point);
        break;
      case State::jarqueBera:  // Compute kurtosis and skewness.
        projection  = dot(point, v) / v_norm;
        sum_first  += projection;
        dummy_term = projection * projection;
        sum_second += dummy_term;
        dummy_term *= projection;
        sum_third  += dummy_term;
        sum_fourth += dummy_term * projection;
        break;
      case State::activeParent:  // Filter the data onwards.
        if (child1->GetDistance(point) < child2->GetDistance(point))
          child1->AddItem(point);
        else
          child2->AddItem(point);
        break;
      // Nothing happens for other two cases.
    }
  }

  // Combines the statistics for two clusters. These clusters should not have
  // two different centers; rather, they should be two equivalent clusters from
  // two different states of the GLA below. Like AddItem, the name is for the
  // sake of consistency and is not relevant to Datapath - this function is only
  // ever called by AddState of the GLA below.
  void AddState(<?=$className?>Cluster & other) {
    switch (state) {
      case State::kMeans: {  // Typical kMeans AddState
        count += other.count;
        sum   += other.sum;
        break;
      }
      case State::gMeans:  // Add to the PVA matrix.
        XTX += other.XTX;
        break;
      case State::kMeansChildren:  // Sum the 2 K-means children.
        child1->AddState(*other.child1);
        child2->AddState(*other.child2);
        break;
      case State::jarqueBera:  // Compute the sum of moments and counts.
        sum_first += other.sum_first;
        sum_second += other.sum_second;
        sum_third  += other.sum_third;
        sum_fourth += other.sum_fourth;
        break;
      case State::activeParent:  // Filter the state down the tree.
        child1->AddState(*other.child1);
        child2->AddState(*other.child2);
        break;
      // Nothing happens for other two cases.
    }
  }

  // This is called during the ShouldIterate step of the GLA.
  void DoIteration(<?=$className?>Node &node) {
    switch (state) {
      case State::kMeans: {  // Move center to centroid
        node.center = sum / count;
        node.count = count;
        break;
      }
      case State::gMeans: {  // Find eigenvalues and find PCA
        // TODO: See if power iteration or armadillo function is faster.
        // TODO: See when DC is faster (probably doesn't matter).
        vec eigenvalues;
        mat eigenvectors;
        uword index;
        eig_sym(eigenvalues, eigenvectors, XTX);
        eigenvalues = abs(eigenvalues);
        double max_eigenvalue = eigenvalues.max(index);
        vec max_eigenvector = eigenvectors.col(index);
        // Norm probably not needed, just a safeguard.
        v = max_eigenvector * sqrt(max_eigenvalue / (count > 1 ? count - 1 : 1))
          / (norm(max_eigenvector, <?=$pMinkowski?>) * ROOT_HALF_PI);
        node.child1 = new <?=$className?>Node(center + v);
        node.child2 = new <?=$className?>Node(center - v);
        node.state = State::kMeansChildren;
        break;
      }
      case State::kMeansChildren: {  // Check for convergence and advance k-means.
        vec previous_center1 = child1->center;
        vec previous_center2 = child2->center;
        child1->DoIteration(*node.child1);
        child2->DoIteration(*node.child2);
        // TODO: Maybe make one child stop if it converges and the other doesn't
        if (   HasConverged(previous_center1, node.child1->center)
            && HasConverged(previous_center2, node.child2->center)) {
          node.state = State::jarqueBera;
          node.child1->state = State::gMeans;
          node.child2->state = State::gMeans;
        } else {
          node.child1->count = 0;
          node.child2->count = 0;
        }
        break;
      }
      case State::jarqueBera: { // Compute kurtosis and skewness; run test.
        // TODO: Maybe use D'Agostino's K-squared Test (or let the user choose).
        long double variance;
        long double kurtosis;
        long double skewness;
        long double JB_statistic;
        long double p_value;
        double mean = sum_first / count;
        variance = sum_second
                 - 2 * mean * sum_first
                 + pow(mean, 2);
        kurtosis = sum_fourth
                 - 4 * mean * sum_third
                 + 6 * pow(mean, 2) * sum_second
                 - 4 * pow(mean, 3) * sum_first
                 + pow(mean, 4);
        skewness = sum_third
                 - 3 * mean * sum_second
                 + 3 * pow(mean, 2) * sum_first
                 - pow(mean, 3);
        kurtosis *= count / pow(variance, 2);
        skewness *= sqrt(count) / pow(variance, 1.5);
        JB_statistic = (pow(skewness, 2) + pow(kurtosis - 3, 2) / 4)
                     * (count / 6);
        // This is equivalent to the Chi squared distribution.
        p_value = boost::math::gamma_p(1.0, JB_statistic / 2);
        if (1 - p_value < <?=$alpha?>) {
          // The data data is not normally distributed, keep the sub clusters.
          node.state = State::activeParent;
        } else {
          // The copy from cluster to node will not copy the children.
          node.state = State::noChildren;
          delete node.child1;
          node.child1 = nullptr;
          delete node.child2;
          node.child2 = nullptr;
        }
        break;
      }
      case State::activeParent: { // Filter the state down the tree.
        child1->DoIteration(*node.child1);
        child2->DoIteration(*node.child2);
        break;
      }
      // Nothing happens for other two cases.
    }
  }

  // Calculation of centroid in order to reposition the center of this cluster
  // for the next iteration.
  vec GetResult() {
    return sum / count;
  }

  // GetDistance from a given point to the center of the cluster. Used to decide
  // which cluster to move a point to in the iterative step. Based on the
  // armadillo implementation of norm which is pre-built to calculate Minkowski
  // distance.
  double GetDistance(vec point) {
    return norm(center - point, <?=$pMinkowski?>);
  }

  // Checks if kMeans has converged. Declared here for correct scoping.
  bool HasConverged(vec old_center, vec new_center) {
    return max(abs((old_center - new_center) / old_center)) < <?=$epsilon?>;
  }
};

class <?=$className?>;

<?  $constantState = lookupResource(
        "statistics::G_Means_Constant_State",
        array(
            'className' => $className,
            'dimension' => $dimension,
        )
    ); ?>

class <?=$className?> {
 typedef <?=$className?>State State;

 private:
  // The typical constant state for an iterable GLA.
  const <?=$constantState?> & constant_state;

  // This is the root cluster; the full tree will be initialized with a single
  // call to the constructor due to the recursive structure of the constructor.
  <?=$className?>Cluster root;

  // A vector whose elements will be individually set to contain the current
  // item's data. The data is packaged as a vector due to the reliance on
  // armadillo. This exists as a member variable as opposed to a local variable
  // in order to avoid repeated memory allocation.
  vec item;

 public:
  <?=$className?>(const <?=$constantState?> & state)
      : constant_state(state),
        root(state.root),
        item(<?=$dimension?>) {
  }

  <?=$className?>(const <?=$className?> &)=delete;
  <?=$className?> & operator =(const <?=$className?> &)=delete;

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
<?  foreach ($codeArray as $name => $counter) { ?>
    item[<?=$counter?>] = <?=$name?>;
<?  } ?>
    // This will automatically make so that the root is set to the center after
    // the first iteration. However, its state needs to be manually changed,
    // which is done in ShouldIterate.
    root.AddItem(item);
  }

  void AddState(<?=$className?> &other) {
    root.AddState(other.root);
  }

  bool ShouldIterate(<?=$constantState?> & modible_state) {
    // All the iterations are processed in the clusters and copied to node.
    root.DoIteration(modible_state.root);
    // Manually increment root's state from computing the mean to G-Means
    if (constant_state.iteration == 0)
      modible_state.root.state = State::gMeans;
    modible_state.root.MakeInactive();
    modible_state.iteration ++;
    return modible_state.root.state != State::inactiveParent
           && modible_state.iteration < 30;
  }

  void GetResult(<?=typed_ref_args($outputs)?>) {
    Json::Value result = constant_state.root.GetOutput(0);
    result["iteration"] = (int) constant_state.iteration;
    cout << result << endl;
  }
};

<?
    return array(
        'kind'             => 'GLA',
        'name'             => $className,
        'system_headers'   => array('armadillo', 'vector', 'string', 'iostream',
                                    'boost/math/special_functions/gamma.hpp',
                                    'boost/math/constants/constants.hpp'),
        'user_headers'     => array(),
        'iterable'         => TRUE,
        'input'            => $inputs,
        'output'           => $outputs,
        'result_type'      => 'single',
        'generated_state'  => $constantState,
        'libraries'        => array('armadillo'),
    );
}
?>
