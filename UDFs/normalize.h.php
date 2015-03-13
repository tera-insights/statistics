<?
// Normalizes an input vector.
function Normalize(array $inputs, array $t_args) {
    grokit_assert(count($inputs) == 1,
                  'Normalize: expected 1 argument.');
    $inputs_ = array_combine(['vector'], array_values($inputs));
    $result = $inputs_['vector'];
    grokit_assert($result->is('vector'),
                  'Normalize: argument 1 should be a vector.');

    $name = generate_name('Normalize');
    $sys_headers  = ['armadillo'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];

?>

<?=$result?> <?=$name?>(<?=const_typed_ref_args($inputs_)?>) {
  return arma::normalise(vector);
}

<?
    return [
        'kind'           => 'FUNCTION',
        'name'           => $name,
        'system_headers' => $sys_headers,
        'user_headers'   => $user_headers,
        'lib_headers'    => $lib_headers,
        'libraries'      => $libraries,
        'input'          => $inputs,
        'result'         => $result,
        'deterministic'  => true,
    ];
}
?>
