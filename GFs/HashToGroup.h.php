<?
function Hash_To_Group_Constant_State(array $t_args) {
    $className = $t_args['className'];
    $states    = $t_args['states'];
    $states_ = array_combine(['hasher'], $states);
    $sys_headers  = [];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = [];
}
?>

class <?=$className?>ConstantState {
 private:
    using Hasher = <?=$states_['hasher']?>;
    const Hasher::Hasher hasher;

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

function Hash_To_Group($t_args, $inputs, $outputs, $states) {
    $className = generate_name('Hash_To_Group');
    $states_ = array_combine(['hasher'], $states);
    $outputs = get_grouping_attributes($t_args);
    $sys_headers  = [];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = [];
    $properties   = ['list'];
    $extra        = [];
    $result_type  = ['multi'];
?>