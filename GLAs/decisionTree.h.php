<?
function Decision_Tree_Constant_State(array $t_args)
{
    // Processing of template arguments.
    $className    = $t_args['className'];
    $domains      = $t_args['domains'];
    $splitType    = get_default($t_args, 'split', 'anova');
    $threshold    = get_default($t_args, 'threshold', 0.01);
    $convLimit    = get_default($t_args, 'conv.limit', 1e-6);
    $emRestarts   = get_default($t_args, 'em.restarts', 3);
    $emDiscIters  = get_default($t_args, 'em.disc.iters', 2);
    $emLearnIters = get_default($t_args, 'em.learn.iters', 30);
    $emMaxTrials  = get_default($t_args, 'em.max.trials', 3);
    $contSplit    = $t_args['cont.split'];
    $contReg      = $t_args['cont.reg'];

    // Processing of local variables.
    $init = implode(', ', $domains);
    $splitMap = ['anova' => 0, 'lda' => 1, 'qda' => 2];
    $split = $splitMap[$splitType];
?>

class <?=$className?>ConstantState {
 public:
	using dist_type = CLUS::Multinormal;
	using splitter_type = CLUS::BinaryObliqueSplitter;
	using tree_type = CLUS::BinaryRegressionTree<dist_type, splitter_type>;
	using size_type = tree_type::SizeType;

 private:
  // The actual tree
	tree_type tree;

  // Iteration counter.
  int iter_counter;

	// Constructor arguments for tree

	// Important: Continuous values vector should contain split variables
	// then regression variables.
	// Size of continuous values vector should be n_cont_split + n_cont_reg
	// # Continuous split variables
	static const constexpr size_type n_cont_split = <?=$contSplit?>;
	// # Continuous regression variables
	static const constexpr size_type n_cont_reg = <?=$contReg?>;

	// 0 = ANOVA, 1 = LDA, 2 = QDA
	static const constexpr int split_type = <?=$split?>;

	// Probability threshold to be considered as part of a node.
	static const constexpr double threshold = <?=$threshold?>;

	// Log likelihood convergence limit for EM
	static const constexpr double conv_limit = <?=$convLimit?>;

	// # EM random restarts
	static const constexpr size_type em_restarts = <?=$emRestarts?>;
	// # Iterations per random restart
	static const constexpr size_type em_disc_iters = <?=$emDiscIters?>;
	// # Iterations during main learning phase
	static const constexpr size_type em_learn_iters = <?=$emLearnIters?>;
	// Max # trials
	static const constexpr size_type em_max_trials = <?=$emMaxTrials?>;

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState()
      : tree(arma::ivec({<?=$init?>}), n_cont_split, n_cont_reg, split_type, threshold,
             conv_limit, em_restarts, em_disc_iters, em_learn_iters, em_max_trials) {
  }
};
<?  return [
        'kind' => 'RESOURCE',
        'name' => $className . 'ConstantState',
        'system_headers' => [],
        'user_headers' => [],
        'lib_headers' => ['regressiontree.h', 'binaryregressiontree.h',
                          'binaryobliquesplitter.h', 'multinormal.h'],
    ];
}

function Decision_Tree(array $t_args, array $inputs, array $outputs)
{
    // Setting output type
    array_set_index($outputs, 0, lookupType('JSON'));

    // Class name is randomly generated.
    $className = generate_name("DTree");

    // Creates an array listing the cardinality of each discrete input.
    $cardinality = function($input) { return $input->get('cardinality'); };
    $domains = array_map($cardinality, array_get_index($inputs, 0)->get('inputs'));

    $sys_headers = ['armadillo', 'limits'];
    $user_headers = [];
    $lib_headers = ['regressiontree.h', 'binaryregressiontree.h',
                    'binaryobliquesplitter.h', 'multinormal.h'];
    $libraries = ['armadillo', 'jsoncpp'];
?>

class <?=$className?>;

<?  $constantState = lookupResource(
        "statistics::Decision_Tree_Constant_State",
        array_merge($t_args,
                    ['className' => $className, 'domains' => $domains]));?>

class <?=$className?> {
 private:
	using tree_type = <?=$constantState?>::tree_type;
	// Type returned by MakeVector on discrete variables
	using dvec = <?=array_get_index($inputs, 0)?>;
	// Type returned by MakeVector on continuous variables
	using cvec = <?=array_get_index($inputs, 1)?>;

	tree_type tree;

 public:
	<?=$className?>(const <?=$constantState?>& const_state)
      : tree(const_state.tree) {
		// GLA States are constructed during a new round, so
		// let the tree know about the new round.
		tree.StartLearningRound();
	}

  // Response variable is appended to the continuous predictors.
	void AddItem(<?=const_typed_ref_args($inputs)?>) {
		tree.LearnSample(<?=array_keys($inputs)[0]?>, <?=array_keys($inputs)[1]?>);
	}

	void AddState(const <?=$className?>& other) {
		tree.Merge(other.tree);
	}

	bool ShouldIterate(<?=$constantState?>& state) {
		// StopLearningRound returns true if the tree is done learning.
    bool should_iterate = tree.StopLearningRound();
    state.tree = tree;
    state.iter_counter++;
    std::cout << "End of iteration " << state.iter_counter << " completed." << std::endl << std::endl;
    int nodes, leaves;
    nodes = leaves = 0;
    state.tree.ComputeSizesTree(nodes, leaves);
    // std::cout << "Nodes: " << nodes << " Leaves: " << leaves << std::endl;
		return !should_iterate;
	}

	void GetResult(<?=typed_ref_args($outputs)?>) {
    Json::Value temp = tree.ToJson();
		<?=array_keys($outputs)[0]?>.set(temp);
	}
};

<?  return [
        'kind'            => 'GLA',
        'name'            => $className,
        'system_headers'  => $sys_headers,
        'user_headers'    => $user_headers,
        'lib_headers'     => $lib_headers,
        'libraries'       => $libraries,
        'iterable'        => true,
        'input'           => $inputs,
        'output'          => $outputs,
        'generated_state' => $constantState,
        'result_type'     => 'single',
    ];
} ?>
