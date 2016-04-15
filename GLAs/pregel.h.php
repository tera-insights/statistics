<?
// This GLA emulates the Pregel operator. During each iteration, each edge in
// the graph is streamed through. The GLA takes in a user-defined function that
// describes how to processes each edge and which messages to send to which
// nodes. The GLA concludes when no messages are passed during an iteration.

// The state passed into the GLA should be one created by the Pregel_Setup GLA.
// It describes the initial value for each vertex, which the messages will then
// update. See said GLA for more information.

// The output is the vertex set with various attached properties.

// Resources:
function Pregel($t_args, $inputs, $outputs, $states)
{
    // Class name is randomly generated.
    $className = generate_name('Pregel');

    // Initializiation of argument names.
    $inputs_ = array_combine(['s', 't', 'weight', 'init'], $inputs);
    $vertex = $inputs_['s'];

    // Initialization of local variables from template arguments.
    $debug = get_default($t_args, 'debug', 1);

    // Construction of outputs.
    $outputs_ = ['node' => $vertex, 'dist' => lookupType('base::double')];
    $outputs = array_combine(array_keys($outputs), $outputs_);

    $sys_headers  = [];
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
        'statistics::Shortest_Path_Constant_State',
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

  // The value of infinity for a double, the initial distance for each vertex.
  static const constexpr double kInf = std::numeric_limits<double>::infinity();

 private:
  // The typical constant state for an iterable GLA.
  const ConstantState& constant_state;

  // This is a reference to the vertex information. Information from the const
  // state is cast to non-const and then updated in parallel by the various
  // workers for this GLA.
  ConstantState::Vertices& vertices;

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
        vertices(const_cast<decltype(vertices)>(state.vertices)),
        iteration(state.iteration),
        finished(true) {
  }

  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    if (iteration == 0) {
      num_nodes = max((long) max(s, t), num_nodes);
      return;
    } else if (iteration == 1) {
      distance(s) = init;
      finished = false;
    } else {
      double new_dist = distance(s) + weight;
      if (new_dist < distance(t)) {
        finished = false;
        distance(t) = new_dist;
      }
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
      distance.set_size(num_nodes);
      distance.fill(kInf);
      return true;
    } else {
<?  if ($debug > 0) { ?>
      cout << "distance: " << distance.t() << endl;
<?  } ?>
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
    // The ordering of operations is important. Don't change it.
    long first = fragment * (num_nodes / kBlock) / num_fragments * kBlock;
    long final = (fragment == num_fragments - 1)
              ? count - 1
              : (fragment + 1) * (count / kBlock) / num_fragments * kBlock - 1;
    return new Iterator(first, final);
  }

  bool GetNextResult(Iterator* it, <?=typed_ref_args($outputs_)?>) {
    if (it->first > it->second)
      return false;
    node = it->first;
    dist = distance(it->first);
    it->first++;
    return true;
  }
};

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

function Pregel_Constant_State($t_args) {
    // Grabbing variables from $t_args
    $className = $t_args['className'];
    $input = $t_args['state'];
    $states = ['state' => $input];
?>

class <?=$className?>ConstantState {
 public:
  using InputState = <?=$input?>;
  using InputTuple = InputState::Tuple;
  using Vertices = InputState::Vector;

 private:
  // The current iteration.
  int iteration;

  // The vertex information. Only one copy is needed.
  const Vertices& vertices;

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState(<?=const_typed_ref_args($state)?>)
      : iteration(0),
        vertices(state.GetVertices()) {
  }
};

<?
    return [
        'kind' => 'RESOURCE',
        'name' => $className . 'ConstantState',
    ];
}

// This GLA is only intended to be used as the input state to the Pregel GLA. It
// simply gathers various attributes associated with a graph vertex and stores
// them in a vector using the vertex ID as the key.

// Each vertex should only be seen once. The vertices are processed in parallel
// but the attributes are updated into a single data structure without locks.

// The output is simply the input.

// Resources:
// tuple: tuple, make_tuple
// vector: vector
function Pregel_Setup($t_args, $inputs, $outputs) {
    // Class name is randomly generated.
    $className = generate_name('PregelState');

    // Initialization of local variables from input names.
    foreach (array_keys(array_slice($inputs, 1)) as $index => $key)
        $vals["val_$index"] = $inputs[$key];
    $inputs_ = array_combine(array_merge(['key'], array_keys($vals)), $inputs);
    $keyType = $inputs_['key'];

    // Processing the outputs.
    $outputs_ = $inputs_;
    $outputs = array_combine(array_keys($outputs), $outputs_);

    $sys_headers  = ['vector', 'tuple'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = [];
    $properties   = [];
    $extras       = [];
?>

class <?=$className?>;

<?  $constantState = lookupResource(
        'statistics::Pregel_Setup_Constant_State',
        ['className' => $className]
    ); ?>

class <?=$className?> {
 public:
  using Tuple = std::tuple<<?=typed($vals)?>>;

  using Vector = std::vector<Tuple>;

 private:
  // The constant state for this GLA.
  using ConstantState = <?=$constantState?>;

  // The gathered information for the vertices.
  static Vector vertices;

  // The number of unique nodes seen.
  long num_nodes;

  // The current iteration.
  int iteration;

  // Used for iterating over the container during GetNextResult.
  int return_counter;

public:
  <?=$className?>()
      : iteration(0),
        num_nodes(state.num_nodes),
        iteration(state.iteration) {
  }

  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    if (iteration == 0)
      num_nodes = max((long) key, num_nodes);
    else
      vertices[key] = std::make_tuple(<?=args($vals)?>);
  }

  void AddState(<?=$className?>& other) {
    if (iteration == 0)
      num_nodes = max(num_nodes, other.num_nodes);
  }

  bool ShouldIterate(ConstantState& state) {
    state.iteration = ++iteration;
<?  if ($debug > 0) { ?>
    cout << "Finished iteration " << iteration << endl;
    cout << "Number nodes: " << num_nodes << endl;
<?  } ?>
    if (iteration == 1) {
      // num_nodes is incremented because IDs are 0-based.
      state.num_nodes = ++num_nodes;
      vertices.reserve(num_nodes);
      return true;
    } else {
      return false;
    }
  }

  void Finalize() {
    return_counter = 0;
  }

  bool GetResult(<?=typed_ref_args($outputs_)?>) const {
    if (return_counter >= vertices.size())
      return false;
    key = return_counter;
<?  foreach (array_keys($outputs) as $index => $name) { ?>
    <?=$name?> = std::get<<?=$index?>>(items[return_counter]);
<?  } ?>
    return true;
  }

  int GetCount() const {
    return items.size();
  }

  const Vector& GetVertices() const {
    return items;
  }
};

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
        'result_type'     => 'multi',
        'generated_state' => $constantState,
    ];
}

function Pregel_Setup_Constant_State(array $t_args)
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
?>
