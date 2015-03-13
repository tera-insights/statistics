<?
// Copyright 2014 Tera Insights, LLC. All Rights Reserved.

// Constructs a tuple containing the given types.

function MakeTuple($inputs, $t_args) {
    $count = count($inputs);

    grokit_warning_assert($count, 'MakeVector: 0 inputs received.');

    foreach (range(0, $count - 1) as $counter)
      $arguments[] = 'arg' . $counter;

    $inputs_ = array_combine($arguments, $inputs);

    $type = lookupType('statistics::tuple', ['types' => array_values($inputs)]);
    $name = generate_name('MakeTuple');

    $sys_headers  = ['tuple'];
    $user_headers = [];
    $lib_headers  = [];
?>

<?=$type?> <?=$name?>(<?=const_typed_ref_args($inputs_)?>) {
  return <?=$type?>(<?=args($inputs_)?>);
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

