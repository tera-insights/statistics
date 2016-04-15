<?
function Connected_Components_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className = $t_args['className'];
?>

class <?=$className?>ConstantState {
 private:
  // The current iteration.
  int iteration;

  // The number of distinct nodes in the graph.
  long num_nodes;

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState()
      : iteration(0),
        num_nodes(0) {
  }
};

<?
    return [
        'kind' => 'RESOURCE',
        'name' => $className . 'ConstantState',
    ];
}

// This GLA computes the connected components of a graph given the edge set as
// inputs. It uses O(V) space and the worst-case run time is O(D * E), where V
// is the total number of vertices; E, the total number of edges; and D, the
// diameter of the widest component.
// The input should be two integers specifying source and target vertices.
// The output is vertex IDs and their component number, which is the smallest ID
// of all the vertices in its connected component.

// Resources:
// armadillo: various data structures
// algorithm: min
function Connected_Components($t_args, $inputs, $outputs)
{
    // Class name is randomly generated.
    $className = generate_name('ConnectedComponents');

    // Initializiation of argument names.
    $inputs_ = array_combine(['s', 't'], $inputs);
    $vertex = $inputs_['s'];

    // Initialization of local variables from template arguments.
    $debug = get_default($t_args, 'debug', 1);

    // Construction of outputs.
    $outputs_ = ['node' => $vertex, 'comp' => $vertex];
    $outputs = array_combine(array_keys($outputs), $outputs_);

    $sys_headers  = ['armadillo', 'algorithm'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $properties   = [];
    $extra        = [];
    $result_type  = ['fragment'];
?>

using namespace arma;
using namespace std;

class <?=$className?>;

<?  $constantState = lookupResource(
        'statistics::Connected_Components_Constant_State',
        ['className' => $className]
    ); ?>

class <?=$className?> {
 public:
  // The constant state for this GLA.
  using ConstantState = <?=$constantState?>;

  // The current and final indices of the result for the given fragment.
  using Iterator = std::pair<int, int>;

  // The work is split into chunks of this size before being partitioned.
  static const constexpr int kBlock = 32;

  // The maximum number of fragments to use.
  static const constexpr int kMaxFragments = 64;

 private:
  // The component for each vertex.
  static arma::Col<<?=$vertex?>> component;

  // The typical constant state for an iterable GLA.
  const ConstantState& constant_state;

  // The number of unique nodes seen.
  long num_nodes;

  // The current iteration.
  int iteration;

  // The number of fragmetns for the result.
  int num_fragments;

  // Whether the algorithm has concluded.
  bool finished;

 public:
  <?=$className?>(const <?=$constantState?>& state)
      : constant_state(state),
        num_nodes(state.num_nodes),
        iteration(state.iteration),
        finished(true) {
  }

  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    if (iteration == 0) {
      num_nodes = max((long) max(s, t), num_nodes);
      return;
    } else {
      // finished flips to false if any component is updated.
      if (finished && component(s) != component(t))
        finished = false;
      // Update the larger component with the smaller.
      component(s) = component(t) = min(component(s), component(t));
    }
  }

  void AddState(<?=$className?> &other) {
    if (iteration == 0)
      num_nodes = max(num_nodes, other.num_nodes);
    else
      finished = finished && other.finished;
  }

  bool ShouldIterate(ConstantState& state) {
    state.iteration = ++iteration;
<?  if ($debug > 0) { ?>
    cout << "finished iteration " << iteration << endl;
<?  } ?>
    if (iteration == 1) {
      // num_nodes is incremented because IDs are 0-based.
      state.num_nodes = ++num_nodes;
<?  if ($debug > 0) { ?>
      cout << "num_nodes: " << num_nodes << endl;
<?  } ?>
      // Allocating space can't be parallelized.
      component = arma::regspace<decltype(component)>(0, num_nodes - 1);
      return true;
    } else {
      return !finished;
    }
  }

  int GetNumFragments() {
    long size = (num_nodes - 1) / kBlock + 1;  // num_nodes / kBlock rounded up.
    num_fragments = (iteration == 0) ? 0 : min(size, (long) kMaxFragments);
<?  if ($debug > 0) { ?>
    cout << "# nodes: " << num_nodes << endl;
    cout << "kBlock: " << kBlock << endl;
    cout << "# fragments: " << size << endl;
    cout << "Returning " << num_fragments << " fragments" << endl;
<?  } ?>
    return num_fragments;
  }

  // The fragment is given as a long so that the below formulas don't overflow.
  Iterator* Finalize(long fragment) {
    long count = num_nodes;
    if (fragment == 1)
      cout << "Count: " << count << endl;
    // The ordering of operations is important. Don't change it.
    long first = fragment * (count / kBlock) / num_fragments * kBlock;
    long final = (fragment == num_fragments - 1)
               ? count - 1
               : (fragment + 1) * (count / kBlock) / num_fragments * kBlock - 1;
<?  if ($debug > 0) { ?>
    printf("Fragment %ld: %ld - %ld\n", fragment, first, final);
<?  } ?>
    return new Iterator(first, final);
  }

  bool GetNextResult(Iterator* it, <?=typed_ref_args($outputs_)?>) {
    if (!finished)
      return false;
    if (it->first > it->second)
      return false;
    node = it->first;
    if (it->first >= component.n_elem)
      cout << "Illegal access. " << it->first << " / " << component.n_elem << endl;
    comp = component(it->first);
    it->first++;
    return true;
  }
};

// Initialize the static member types.
arma::Col<<?=$vertex?>> <?=$className?>::component;

typedef <?=$className?>::Iterator <?=$className?>_Iterator;

<?
    return [
        'kind'            => 'GLA',
        'name'            => $className,
        'system_headers'  => $sys_headers,
        'user_headers'    => $user_headers,
        'lib_headers'     => $lib_headers,
        'libraries'       => $libraries,
        'properties'      => $properties,
        'extra'           => $extra,
        'iterable'        => true,
        'intermediates'   => true,
        'input'           => $inputs,
        'output'          => $outputs,
        'result_type'     => $result_type,
        'generated_state' => $constantState,
    ];
}
?>
