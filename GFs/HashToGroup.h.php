<?
function Hash_To_Group_Constant_State(array $t_args) {
    $className = $t_args['className'];
    $states    = $t_args['states'];
    $states_ = array_combine(['hasher'], $states);
    $sys_headers  = [];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = [];
?>

class <?=$className?>ConstantState {
 private:
    using Hasher = <?=$states_['hasher']?>;
    const Hasher hasher;

 public:
    friend class <?=$className?>;

    <?=$className?>ConstantState(<?=const_typed_ref_args($states_)?>)
      : hasher(hasher) {
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

function get_grouping_attributes_for_hash_to_group(array $t_args, array $inputs) {
    $groupingInputNames = $t_args['group'];
    $grouping_attributes = [];
    foreach( $inputs as $name => $type ) {
        $is_grouping_attribute = in_array($name, $groupingInputNames);
        if($is_grouping_attribute) {
        $grouping_attributes[$name] = $type;
        }
    }
    return $grouping_attributes;
}


function Hash_To_Group(array $t_args, array $inputs, array $states) {
    $className = generate_name('Hash_To_Group');
    $states_ = array_combine(['hasher'], $states);
    $groupingInputs = get_grouping_attributes_for_hash_to_group($t_args, $inputs);
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
        "statistics::Hash_To_Group_Constant_State", ['className' => $className, 'states' => $states]
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
    return !constant_state.hasher.IsGroupSurvivor(<?=args($groupingInputs)?>);
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


