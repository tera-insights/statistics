<?
// Copyright 2014 Tera Insights, LLC. All Rights Reserved.

// Constructs a tuple containing the given types.

function MakeVector($inputs, $t_args) {
    grokit_assert(count($inputs) > 0, 'MakeVector: 0 inputs received.');

    foreach (range(0, $count - 1) as $counter)
      $arguments[] = 'arg' . $counter;

    $inputs_ = array_combine($arguments, $inputs);

    $type = lookupType('statistics::tuple', ['types' => array_values($inputs)]);
    $name = generate_name('MakeTuple');

    $sys_headers  = ['tuple'];
    $user_headers = [];
    $lib_headers  = [];
?>

<?=$type?> <?=$funName?>(<?=const_typed_ref_args($inputs_)?>) {
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
        'result'         => $output,
        'deterministic'  => true,
    ];
}
?>

