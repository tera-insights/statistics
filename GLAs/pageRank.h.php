<?
function Page_Rank_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className          = $t_args['className'];
?>
using namespace arma;
using namespace std;

class <?=$className?>ConstantState {
 private:
  // The current iteration.
  int iteration;

  // The number of distinct nodes in the graph.
  int num_nodes;

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

// This GLA computes the page ranks of a graph given the edge set as inputs. The
// algorithm uses O(V) time and takes O(I * E) time, where V is the total number
// of vertices; E, the total number of edges; and I, the number of iterations.
// The input should be two integers specifying source and target vertices.
// The output is vertex IDs and their page rank, expressed as a float.
// Template Args:
// adj: Whether the edge count and page rank of a given vertex should be stored
//   adjacently. Doing so reduced random lookups but increases update time.
// hash: Whether the key IDs need to be converted to zero-based indices.
// Resources:
// armadillo, vector, unordered_map: various data structures.
function Page_Rank($t_args, $inputs, $outputs)
{
    // Class name is randomly generated.
    $className = generate_name('PageRank');

    // Initializiation of argument names.
    $inputs_ = array_combine(['s', 't'], $inputs);
    $vertex = $inputs_['s'];

    // Initialization of local variables from template arguments.
    $adj = $t_args['adj'];

    // Construction of outputs.
    $outputs_ = ['node' => $vertex, 'rank' => lookupType('float')];
    $outputs = array_combine(array_keys($outputs), $outputs_);

    $sys_headers  = ['armadillo', 'vector', 'mct/hash-map.hpp', 'unordered_map', 'algorithm'];
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
        'statistics::Page_Rank_Constant_State',
        ['className' => $className]
    ); ?>

class <?=$className?> {
 public:
  // The constant state for this GLA.
  using ConstantState = <?=$constantState?>;

  // The key type for the map.
  using Key = uint64_t;

  // A mapping of vertex IDs to a zero-based enumeration.
  using Map = mct::closed_hash_map<Key, uint64_t>;

  // The list of keys, a reverse mapping of the above.
  using KeySet = std::vector<Key>;

  // The current and final indices of the result for the given fragment.
  using Iterator = std::pair<int, int>;

  // The type of the vertex IDs.
  using Vertex = <?=$vertex?>;

  // The factor by which to scale the size of the vectors;
  static const constexpr int kScale = 2;

  // The size various objects are initialized to have.
  static const constexpr int kSize = 100;

  // The value of the damping constant used in the page rank algorithm.
  static const constexpr float kDamping = 0.85;

  // The number of iterations to perform, not counting the initial set-up.
  static const constexpr int kIterations = 20;

  // The work is split into chunks of this size before being partitioned.
  static const constexpr int kBlock = 32;

  // The maximum number of fragments to use.
  static const constexpr int kMaxFragments = 64;

 private:
<?  if ($adj) { ?>
  // The edge weighting and page rank are stored side-by-side. This will make
  // updating the page ranks more costly but reduces random lookups. The weight
  // is the precomputed inverse of the number of edges for that vertex.
  static arma::mat info;
<?  } else { ?>
  // The edge weighting for each vertex. The weight is the precomputed inverse
  // of the number of edges for that vertex.
  static arma::rowvec weight;

  // The page-rank for each vertex.
  static arma::rowvec rank;
<?  } ?>

  // The value of the summation over adjacent nodes for each vertex.
  static arma::rowvec sum;

  // The set of vertices are mapped to zero-based consecutive indices.
  static Map index_map;

  // A reverse mapping of the above, in which the key is just the index.
  static KeySet key_set;

  // The typical constant state for an iterable GLA.
  const ConstantState& constant_state;

  // The local objects used to compute the above mappings.
  Map indices;
  KeySet keys;

  // The number of unique nodes seen.
  long num_nodes;

  // The current iteration.
  int iteration;

  // The number of fragmetns for the result.
  int num_fragments;

 public:
  <?=$className?>(const <?=$constantState?>& state)
      : constant_state(state),
        num_nodes(state.num_nodes),
        iteration(state.iteration) {
  }

  // Basic dynamic array allocation.
  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
    if (iteration == 0) {
      Update(s, true);
      Update(t, false);
      return;
    }
    uint64_t s_index = index_map[Hash(s)];
    if (iteration == 1) {
<?  if ($adj) { ?>
      info(1, s_index)++;
<?  } else { ?>
      weight(s_index)++;
<?  } ?>
      return;
    }
    uint64_t t_index = index_map[Hash(t)];
<?  if ($adj) { ?>
    sum(t_index) += prod(info.col(s_index));
<?  } else { ?>
    sum(t_index) += weight(s_index) * rank(s_index);
<?  } ?>
  }

  // Hashes are merged.
  void AddState(<?=$className?> &other) {
    if (iteration == 0) {
      for (auto it = other.indices.begin(); it != other.indices.end(); ++it) {
        // Iterate over the vertices seen by the other state.
        const Key& key = it->first;
        auto match = indices.find(key);
        if (match == indices.end()) {
          // This state has not seen the current vertex.
          indices.insert(make_pair(key, num_nodes++));
          keys.push_back(key);
        }
      }
    }
  }

  // Most computation that happens at the end of each iteration is parallelized
  // by performed it inside Finalize.
  bool ShouldIterate(ConstantState& state) {
    state.iteration = ++iteration;
    cout << "finished iteration " << iteration << endl;
    if (iteration == 1) {
      state.num_nodes = num_nodes;
      cout << "num_nodes: " << num_nodes << endl;
      // These operation are not easily parallelized and are performed here.
      index_map.swap(indices);
      key_set.swap(keys);
      // Allocating space can't be parallelized.
      sum.set_size(num_nodes);
<?  if ($adj) { ?>
      info.set_size(2, num_nodes);
      info.row(0).fill(1);
<?  } else { ?>
      rank.set_size(num_nodes);
      weight.set_size(num_nodes);
      rank.fill(1);
<?  } ?>
      return true;
    } else if (iteration == 2) {
<?  if ($adj) { ?>
      info.row(1) = pow(info.row(1), -1);
<?  } else { ?>
      weight = pow(weight, -1);
<?  } ?>
      return true;
    } else {
      cout << "sum: " << accu(sum) << endl;
      cout << "pr: " << accu(info.row(0)) << endl;
      return iteration < kIterations + 1;
    }
  }

  int GetNumFragments() {
    long size = (num_nodes - 1) / kBlock + 1;  // num_nodes / kBlock rounded up.
    num_fragments = (iteration <= 1) ? 0 : min(size, (long) kMaxFragments);
    cout << "Returning " << num_fragments << " fragments" << endl;
    return num_fragments;
  }

  Iterator* Finalize(int fragment) {
    int count = key_set.size();
    // The ordering of operations is important. Don't change it.
    int first = fragment * (count / kBlock) / num_fragments * kBlock;
    int final = (fragment == num_fragments - 1)
              ? count - 1
              : (fragment + 1) * (count / kBlock) / num_fragments * kBlock - 1;
    // printf("fragment: %d\tfirst: %d\tfinal: %d\n", fragment, first, final);
<?  if ($adj) { ?>
    info.row(0).subvec(first, final) = (1 - kDamping)
                                     + kDamping * sum.subvec(first, final);
<?  } else { ?>
    rank.subvec(first, final) = (1 - kDamping)
                              + kDamping * sum.subvec(first, final);
<?  } ?>
    sum.subvec(first, final).zeros();
    return new Iterator(first, final);
  }

  bool GetNextResult(Iterator* it, <?=typed_ref_args($outputs_)?>) {
    if (iteration < kIterations + 1)
      return false;
    if (it->first > it->second)
      return false;
    node = key_set[it->first];
<?  if ($adj) { ?>
    rank = info(0, it->first);
<?  } else { ?>
    rank = rank(it->first);
<?  } ?>
    it->first++;
    return true;
  }

 private:
  // Updates the graph information based on vertex and if the edge it was part
  // of originated from it.
  void Update(Vertex v, bool outgoing) {
    Key key = Hash(v);
    auto it = indices.find(key);
    if (it == indices.end()) {
      // New node seen. Add it to map.
      indices.insert(make_pair(key, num_nodes++));
      keys.push_back(key);
    }
  }
};

// Initialize the static member types.
<?  if ($adj) { ?>
arma::mat <?=$className?>::info;
<?  } else { ?>
arma::rowvec <?=$className?>::weight;
arma::rowvec <?=$className?>::rank;
<?  } ?>
arma::rowvec <?=$className?>::sum;
<?=$className?>::Map <?=$className?>::index_map;
<?=$className?>::KeySet <?=$className?>::key_set;


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
