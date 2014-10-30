<?
// Normalizes an input vector.
function Normalize(array $inputs, array $t_args) {
    grokit_assert(count($inputs) == 1,
                  "Normalize must be given a single argument.");
    $name = "input";
    $inputs = array_combine(["input"], array_values($inputs));
    $input = $inputs[$name];
    grokit_assert($input->is("array"), "Normalize must be given an array type.");
    $size = $input->get("size");
    $type = $input->get("type");
    grokit_assert($type->is("numeric"), "Only numeric arrays can be normalized.");

    $funName = generate_name("Normalize");
    $sys_headers = ['cmath'];
    $user_headers = [];
    $lib_headers = [];
    $output = lookupType("base::fixedarray", ['size' => $size, 'type' => lookupType("base::double")]);

?>

<?=$output?> <?=$funName?>(<?=const_typed_ref_args($inputs)?>) {
  using namespace std;
  using namespace arma;

  double l2_norm = 0.0;
  for (<?=$type?> element : <?=$name?>)
    l2_norm += element * element;
  l2_norm = sqrt(l2_norm);
  <?=$output?> value;
  for (int counter = 0; counter < <?=$size?>; counter++)
    value[counter] = <?=$name?>[counter] / l2_norm;
  return value;

  // Col<<?=$type?>> copy(<?=$name?>.data(), <?=$size?>);
  // <?=$output?> value;
  // vec::fixed<<?=$size?>>  normalized = normalise(copy);
  // value.from_memory(normalized.memptr());
  // return value;
}

<?
    return [
        'kind'           => 'FUNCTION',
        'name'           => $funName,
        'system_headers' => $sys_headers,
        'user_headers'   => $user_headers,
        'lib_headers'    => $lib_headers,
        'input'          => $inputs,
        'result'         => $output,
        'deterministic'  => true,
    ];
}
?>
