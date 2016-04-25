<?
// This GLA emulates the Pregel operator. During each iteration, each edge in
// the graph is streamed through. The GLA takes in a user-defined function that
// describes how to processes each edge and which messages to send to which
// nodes. The GLA concludes when no messages are passed during an iteration.

// The state passed into the GLA should be one created by the Pregel_Setup GLA.
// It describes the initial value for each vertex, which the messages will then
// update. See said GLA for more information.

// The output is the vertex set with various attached properties.

// Template Args:
// properties: The names of the edge properties.
// vertices: The name of the two vertices for an edge.
// combine: The specification for how to update a vertex' properties given
//   another vertex and their edge's properties. It must return a boolean
//   describing whether an update was performed.
// directed: Whether to only update the target vertex.

// Resources:
function Pregel($t_args, $inputs, $outputs, $states)
{
    // Class name is randomly generated.
    $className = generate_name('Pregel');

    // Initialization of local variables from template arguments.
    $directed = get_default($t_args, 'directed', true);
    $debug = get_default($t_args, 'debug', 1);
    $prop = $t_args['properties'];
    grokit_assert(is_array($prop), 'attributes should be a string array');
    $numProp = count($prop);
    $expProp = count($inputs) - 2;
    grokit_assert($numProp == $expProp,
                  "Expected $expProp attribututes. Got $numProp.");
    $vertices = array_values($t_args['vertices']);
    $vertexName = $t_args['name'];
    $sep = $numProp ? ', ' : '';  // The separator before the property args.
    $message = $t_args['message'];
    $combine = $t_args['combine'];
    $intermediates = $combine != "";

    // Initializiation of argument names.
    $inputs_ = array_combine(array_merge(['s', 't'], $prop), $inputs);
    $prop = array_combine($prop, array_slice($inputs, 2));
    $vertex = $inputs_['s'];

    // Processing of input state.
    $states_ = array_combine(['state'], $states);
    $state = $states_['state'];
    $atts = $state->get('attributes');

    // Construction of outputs.
    $outputs_ = array_merge(['vertex' => $vertex], $atts);
    $outputs = array_combine(array_keys($outputs), $outputs_);

    $sys_headers  = ['armadillo'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $properties   = [];
    $extra        = [];
    $result_type  = ['fragment'];
?>

using namespace std;

class <?=$className?>;

<?  $constantState = lookupResource(
        'statistics::Pregel_Constant_State',
        ['className' => $className, 'state' => $state]
    ); ?>

class <?=$className?> {
 public:
  // The constant state for this GLA.
  using ConstantState = <?=$constantState?>;

  // The type for each vertex.
  using Vertex = ConstantState::InputVertex;

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
        num_nodes(vertices.size()),
        iteration(state.iteration),
        finished(true) {
    auto state_copy = const_cast<<?=$constantState?>&>(state);
    cout << "Time taken for last iteration: " << state_copy.timer.toc() << endl;
  }

  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    // printf("Processing edge: %u -> %u\n", s, t);
    finished = finished & Message(vertices[t], vertices[s]<?=$sep, args($prop)?>);
<?  if (!$directed) { ?>
    finished = finished & Message(vertices[s], vertices[t]<?=$sep, args($prop)?>);
<?  } ?>
    // if (finished) cout << "Update performed." << endl;
  }

  void AddState(<?=$className?> &other) {
    finished = finished && other.finished;
  }

  bool ShouldIterate(ConstantState& state) {
    state.iteration = ++iteration;
<?  if ($debug > 0) { ?>
    cout << "finished iteration " << iteration << endl;
<?  } ?>
    return !finished;
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
    vertex = it->first;
<?  if ($intermediates) { ?>
    Combine(vertices[vertex]);
    if (!finished)
      return false;
<?  } ?>
<?  foreach (array_keys($atts) as $name) { ?>
    <?=$name?> = vertices[vertex].<?=$name?>;
<?  } ?>
    it->first++;
    return true;
  }

 private:
  // This updates the target vertex given the source vertex and the properties
  // of their shared edge. It returns whether an update was actually performed.
  bool Message(Vertex& <?=$vertices[1]?>, Vertex& <?=$vertices[0]?>
              <?=$sep, const_typed_ref_args($prop)?>) {
    // cout << "Old weight: " << T.D << " New weight: " << (S.D + W) << endl;
    <?=str_replace('\n', "\n", $message)?>
  }

  // This is called at the end of each iteration per vertex. It processes the
  // various messages received during an iteration and sets up the vertex for
  // the next iteration.
  void Combine(Vertex& <?=$vertexName?>) {
    <?=str_replace('\n', "\n", $combine)?>
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
        'intermediates'   => $intermediates,
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
  using InputVertex = InputState::Vertex;
  using Vertices = InputState::Vector;

 private:
  // The current iteration.
  int iteration;

  // The vertex information. Only one copy is needed.
  const Vertices& vertices;

  arma::wall_clock timer;

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState(<?=const_typed_ref_args($states)?>)
      : iteration(0),
        vertices(state.GetVertices()) {
    cout << "Input has " << state.GetCount() << " vertices." << endl;
    timer.tic();
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

// Template Args:
// attributes: The names of the vertex properties to store.

// Resources:
// tuple: tuple, make_tuple
// vector: vector
// algorithm: max
function Pregel_Setup($t_args, $inputs, $outputs) {
    // Class name is randomly generated.
    $className = generate_name('PregelState');

    // Initialization of local variables from template arguments.
    $atts = $t_args['attributes'];
    $debug = get_default($t_args, 'debug', 1);
    grokit_assert(is_array($atts), 'attributes should be a string array');
    $numAtts = count($atts);
    $expAtts = count($inputs) - 1;
    grokit_assert($numAtts == $expAtts,
                  "Expected $expAtts attribututes. Got $numAtts.");

    // Initialization of local variables from input names.
    $inputs_ = array_combine(array_merge(['key'], $atts), $inputs);
    $atts = array_combine($atts, array_slice($inputs, 1));
    $keyType = $inputs_['key'];

    // Processing the outputs.
    $outputs_ = $inputs_;
    $outputs = array_combine(array_keys($outputs), $outputs_);

    $sys_headers  = ['vector', 'tuple', 'algorithm', 'armadillo'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $properties   = [];
    $extra        = ['attributes' => $atts];
?>

using namespace std;

class <?=$className?>;

<?  $constantState = lookupResource(
        'statistics::Pregel_Setup_Constant_State',
        ['className' => $className, 'attributes' => $atts]
    ); ?>

class <?=$className?> {
 public:
  using Vertex = <?=$constantState?>::Vertex;
  using Vector = <?=$constantState?>::Vector;

 private:
  // The constant state for this GLA.
  using ConstantState = <?=$constantState?>;

  // The number of unique nodes seen.
  long num_nodes;

  // The current iteration.
  int iteration;

  // The gathered information for the vertices.
  Vector& vertices;

  // Used for iterating over the container during GetNextResult.
  int return_counter;

public:
  <?=$className?>(const <?=$constantState?>& state)
      : num_nodes(state.num_nodes),
        iteration(state.iteration),
        vertices(const_cast<decltype(vertices)>(state.vertices)) {
  }

  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    if (iteration == 0)
      num_nodes = max((long) key, num_nodes);
    else
      vertices[key] = Vertex{<?=args($atts)?>};
  }

  void AddState(<?=$className?>& other) {
    if (iteration == 0)
      num_nodes = max(num_nodes, other.num_nodes);
  }

  bool ShouldIterate(ConstantState& state) {
    state.iteration = ++iteration;
<?  if ($debug > 0) { ?>
    cout << "Setup: Finished iteration " << iteration << endl;
<?  } ?>
    if (iteration == 1) {
      // num_nodes is incremented because IDs are 0-based.
      state.num_nodes = ++num_nodes;
<?  if ($debug > 0) { ?>
      cout << "Setup: Number nodes: " << num_nodes << endl;
<?  } ?>
      vertices.resize(num_nodes);
      return true;
    } else {
<?  if ($debug > 0) { ?>
      cout << "Setup: Finished with number nodes: " << num_nodes << endl;
      cout << "Time taken for last setup: " << state.timer.toc() << endl;
<?  } ?>
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
<?  foreach (array_keys($atts) as $name) { ?>
    <?=$name?> = vertices[return_counter].<?=$name?>;
<?  } ?>
    return true;
  }

  long GetCount() const {
    return vertices.size();
  }

  const Vector& GetVertices() const {
    return vertices;
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
    $atts      = $t_args['attributes'];
?>

class <?=$className?>ConstantState {
 public:
  struct Vertex {
    <?=str_replace(',', ';', typed_args($atts))?>;
  };

  using Vector = std::vector<Vertex>;

 private:
  // The current iteration.
  int iteration;

  // The number of distinct nodes in the graph.
  long num_nodes;

  // Only one set is maintained.
  Vector vertices;

  arma::wall_clock timer;

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState()
      : iteration(0),
        num_nodes(0) {
    timer.tic();
  }
};

<?
    return [
        'kind' => 'RESOURCE',
        'name' => $className . 'ConstantState',
    ];
}
?>
