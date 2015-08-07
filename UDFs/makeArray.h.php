<?
// Copyright 2014 Tera Insights, LLC. All Rights Reserved.

// Constructs a fixed array from a fixed vector.

function MakeArray($inputs, $t_args) {
    $count = count($inputs);

    grokit_assert($count == 1, 'MakeArray: exactly 1 input expected.');

    $inputs_ = array_combine(['vector'], $inputs);

    grokit_assert($inputs_['vector']->is('vector'),
                  'MakeArray: vector input expected.');

    $type = lookupType('base::Array',
                       ['type' => $inputs_['vector']->get('type'),
                        'size' => $inputs_['vector']->get('size')]);
    $name = generate_name('MakeArray');

    $sys_headers  = [];
    $user_headers = [];
    $lib_headers  = [];
?>

<?=$type?> <?=$name?>(<?=const_typed_ref_args($inputs_)?>) {
  <?=$type?> result;
  result.from_memory(vector.memptr());
  return result;
}

<?
    return [
        'kind'           => 'FUNCTION',
        'name'           => $name,
        'system_headers' => $sys_headers,
        'user_headers'   => $user_headers,
        'lib_headers'    => $lib_headers,
        'input'          => $inputs,
        'result'         => $type,
        'deterministic'  => true,
    ];
}
?>

