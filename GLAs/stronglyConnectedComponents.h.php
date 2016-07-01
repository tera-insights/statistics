<?
function Strongly_Connected_Components_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className = $t_args['className'];
?>

// This enumeration represents the phase of the strongly connected components algorithm.
enum class SCCPhase { INITIALIZE, TRIMMING, FORWARD, BACKWARD };

class <?=$className?>ConstantState {
 private:
  // The current iteration.
  int iteration;

  // The number of distinct nodes in the graph.
  long num_nodes;

  // The current phase of the algorithm.
  SCCPhase phase;

 public:
   friend class <?=$className?>;

   <?=$className?>ConstantState()
     : iteration(0),
     num_nodes(0),
     phase(SCCPhase::INITIALIZE) { }
};

<?
    return [
        'kind' => 'RESOURCE',
        'name' => $className . 'ConstantState',
    ];
}

// This GLA computes the strongly connected components of a directed graph given the
// edge set as inputs. It uses O(V) space and the worst-case run time is O(D * E),
// where V is the total number of vertices; E, the total number of edges; and D, the
// diameter of the widest component.

// The input should be two integers specifying source and target vertices.
// The output is vertex IDs and their component number, which is the smallest ID
// of all the vertices in its connected component.

// Resources:
// armadillo: various data structures
// algorithm: min

function Strongly_Connected_Components($t_args, $inputs, $outputs)
{
    // Class name is randomly generated.
    $className = generate_name('StronglyConnectedComponents');
    // Initialization of argument names.
    $inputs_ = array_combine(['s', 't'], $inputs);
    $vertex = $inputs_['s'];
    // Processing of template arguments.
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
        'statistics::Strongly_Connected_Components_Constant_State',
        ['className' => $className]
    ); ?>

class <?=$className?> {
public:
  // The constant state for this GLA.
  using ConstantState = <?=$constantState?>;

  // The current and final indices of the result for the given fragment.
  using Iterator = std::pair<int, int>;

  // The work is split into chunks of this size before being partitioned.
  static constexpr int kBlock = 32;

  // The maximum number of fragments to use.
  static constexpr int kMaxFragments = 64;

private:
  // The component ID for each vertex.
  static arma::Col<<?=$vertex?>> component;

  // Whether a vertex is inactive.
  static arma::Row<int> inactive;

  // Whether the in-degree of a vertex is not zero.
  static arma::Row<int> indegree;

  // Whether the out-degree of a vertex is not zero.
  static arma::Row<int> outdegree;

  // The typical constant state for an iterable GLA.
  const ConstantState& constant_state;

  // The number of unique nodes seen.
  long num_nodes;

  // The current iteration.
  int iteration;

  // The current phase of the algorithm.
  SCCPhase phase;

  // The number of fragments for the result.
  int num_fragments;

  // Whether the a phase of the algorithm has concluded.
  bool finished;

  // Whether the algorithm has concluded
  bool terminated;

public:
  <?=$className?>(const <?=$constantState?>& state)
    : constant_state(state),
      num_nodes(state.num_nodes),
      iteration(state.iteration),
      phase(state.phase),
      finished(true),
      terminated(true) {
  }

  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    if (phase == SCCPhase::INITIALIZE) {
      num_nodes = std::max(static_cast<long>(std::max(s, t)), num_nodes);
    } else if (phase == SCCPhase::TRIMMING) {
      // if both s and t are active, update their in-degree/out-degree
      if (!inactive(s) && !inactive(t)) {
        indegree(t) = 1;
        outdegree(s) = 1;
      }
    } else if (phase == SCCPhase::FORWARD) {
      if (!inactive(s) && !inactive(t)) {
        finished = finished & (component(s) == component(t));
        // Update the larger component with the smaller.
        component(s) = component(t) = std::min(component(s), component(t));
      }
    } else if (phase == SCCPhase::BACKWARD) {
      // if a vertex is active and its component id equals to its vertex id, it's the root of the
      // strongly connected component it belongs to
      if (!inactive(s) && component(s) == s) {
        finished = false;
        inactive(s) = 1;
      }
      if (!inactive(t) && component(t) == t) {
        finished = false;
        inactive(t) = 1;
      }
      // once we make the root of a SCC inactive, find out all other vertices in this component using
      // reverse edges (t -> s)
      if (inactive(t) && !inactive(s) && component(t) == component(s)) {
        finished = false;
        inactive(s) = 1;
      }
    }
  }

  void AddState(<?=$className?> &other) {
    if (phase == SCCPhase::INITIALIZE) {
      num_nodes = std::max(num_nodes, other.num_nodes);
    } else if (phase == SCCPhase::FORWARD || phase == SCCPhase::BACKWARD) {
      finished &= other.finished;
    }
  }

  bool ShouldIterate(ConstantState& state) {
    state.iteration = ++iteration;
<?  if ($debug > 0) { ?>
    cout << "finished iteration " << iteration << endl;
<?  } ?>
    if (phase == SCCPhase::INITIALIZE) {
      state.num_nodes = ++num_nodes;
      component = arma::regspace<decltype(component)>(0, num_nodes - 1);
      inactive.zeros(num_nodes);
      indegree.zeros(num_nodes);
      outdegree.zeros(num_nodes);
      phase = SCCPhase::TRIMMING;
      state.phase = phase;
<?  if ($debug > 0) { ?>
      cout << "finished INITIALIZE phase" << endl;
<?  } ?>
      return true;
    } else if (phase == SCCPhase::TRIMMING) {
      for (int i = 0; i < num_nodes; i++) {
        // identify trivial SCCs (vertices with in-degree or out-degree equal to 0)
        if (indegree(i) == 0 || outdegree(i)== 0) {
          inactive(i) = 1;
        }
      }
      // jump to forward traversal phase
      phase = SCCPhase::FORWARD;
      state.phase = phase;
<?  if ($debug > 0) { ?>
      cout << "finished TRIMMING phase" << endl;
<?  } ?>
      return true;
    } else if (phase == SCCPhase::FORWARD) {
      if (finished) {
        // jump to backward traversal phase
        phase = SCCPhase::BACKWARD;
        state.phase = phase;
<?  if ($debug > 0) { ?>
        cout << "finished FORWARD phase" << endl;
<?  } ?>
      }
      return true;
    } else if (phase == SCCPhase::BACKWARD) {
      if (finished) {
        for (int i = 0; i < num_nodes; i++) {
          // reset attributes if a vertex is active
          if (!inactive(i)) {
            component(i) = i;
            indegree(i) = 0;
            outdegree(i) = 0;
            terminated = false;
          }
        }
        if (terminated) {
          return false;
        }
        // if not terminated, jump back to trimming phase and start another super step
        phase = SCCPhase::TRIMMING;
        state.phase = phase;
<?  if ($debug > 0) { ?>
        cout << "finished BACKWARD phase" << endl;
<?  } ?>
      }
      return true;
    }
    return false;
  }

  int GetNumFragments() {
    long size = (num_nodes - 1) / kBlock + 1;  // num_nodes / kBlock rounded up.
    num_fragments = (iteration == 0) ? 0 : min(size, (long) kMaxFragments);
<?  if ($debug > 0) { ?>
   //cout << "# nodes: " << num_nodes << endl;
   //cout << "kBlock: " << kBlock << endl;
   //cout << "# fragments: " << size << endl;
   //cout << "Returning " << num_fragments << " fragments" << endl;
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
    // printf("Fragment %ld: %ld - %ld\n", fragment, first, final);
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
arma::Row<int> <?=$className?>::inactive;
arma::Row<int> <?=$className?>::indegree;
arma::Row<int> <?=$className?>::outdegree;

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
        'intermediates'   => false,
        'input'           => $inputs,
        'output'          => $outputs,
        'result_type'     => $result_type,
        'generated_state' => $constantState,
    ];
}
?>
