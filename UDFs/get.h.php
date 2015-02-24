<?
// Copyright 2014 Tera Insights, LLC. All Rights Reserved.

// Analagous to std::get for std::tuple.

function MakeVector($inputs, $t_args) {
    grokit_assert(count($inputs) == 1, 'Get: input should be a single tuple.');
    $type = array_get_index($inputs, 0)
    grokit_assert($type->is('tuple'), 'Get: input should be a single tuple.');

    grokit_assert(count($t_args) == 1, 'Get: expected one template argument.');
    $index = array_get_index($t_args, 0);
    grokit_assert(is_int($index), 'Get: index should be an integer.');

    $size = count($type->get('size'));
    grokit_assert(0 <= $index && $index < $size, 'Get: index out of bounds.');

    $inputs_ = array_combine(['tuple'], $inputs);

    $type = array_get_index($type->get('types'), $index);
    $name = generate_name('GetTuple');

    $sys_headers  = ['tuple'];
    $user_headers = [];
    $lib_headers  = [];
?>

<?=$type?> <?=$name?>(<?=const_typed_ref_args($inputs_)?>) {
  return get<<?=$index?>>(tuple);
}

<?
    return [
        'kind'           => 'FUNCTION',
        'name'           => $name,
        'system_headers' => $sys_headers,
        'user_headers'   => $user_headers,
        'lib_headers'    => $lib_headers,
        'input'          => $inputs,
        'result'         => $output,
        'deterministic'  => true,
    ];
}
?>

