<?
function Segmented_Hash_To_Group_Constant_State(array $t_args) {
    $className = $t_args['className'];
    $states    = $t_args['states'];
    $states_ = array_combine(['segmenter'], $states);
    $sys_headers  = [];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = [];
?>

class <?=$className?>ConstantState {
 private:
    using Segmenter = <?=$states_['segmenter']?>;
    const Segmenter &segmenter;
    unsigned long minimum_bucket_score;

 public:
    friend class <?=$className?>;

    <?=$className?>ConstantState(<?=const_typed_ref_args($states_)?>)
      : segmenter(segmenter) {
          minimum_bucket_score = 0;
    }

    const Segmenter &GetSegmenter() {
        return segmenter;
    }
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

function get_segmenting_inputs(array $t_args, array $inputs) {
    return $t_args['segmenting_inputs'];
}

function Segmented_Hash_To_Group(array $t_args, array $inputs, array $states) {
    $className = generate_name('Segmented_Hash_To_Group');
    $states_ = array_combine(['hasher'], $states);
    $groupingInputs = get_grouping_attributes_for_hash_to_group($t_args, $inputs);
    $segmentingInputs = get_segmenting_inputs($t_args, $inputs);
    $sys_headers  = [];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = [];
    $properties   = ['list'];
    $extra        = [];
    $result_type  = ['multi'];
?>

class <?=$className?>;

<?  $constantState = lookupResource(
        "statistics::Segmented_Hash_To_Group_Constant_State", [
            'className' => $className,
            'states' => $states
        ]
    ); ?>

class <?=$className?> {
 private:
  using ConstantState = <?=$constantState?>;
  const <?=$constantState?>& constant_state;

 public:
  <?=$className?>(const <?=$constantState?>& state)
      : constant_state(state) {
  }

  bool Filter(<?=const_typed_ref_args($inputs)?>) {
    auto segment =
        constant_state.segmenter.GetInnerGLA(<?=args($segmentingInputs)?>);
    auto bucket_score = segment.GetBucketScore(<?=args($groupingInputs)?>);
    return bucket_score >= constant_state.minimum_bucket_score;
  }
};

<?
    return [
		'kind'				=> 'GF',
		'name'				=> $className,
		'input'				=> $inputs,
		'system_headers' 	=> $sys_headers,
		'libraries'			=> $libraries,
		'generated_state'	=> $constantState,
	];
}
?>


