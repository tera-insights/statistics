<?
// This GLA builds a fixed-sized multi-dimensional array. Given a list of scalar
// indices followed by a scalar or an array, it accumulates the data into an
// array of higher dimension.

// Aggregate: How to aggregate data when a single point in the result is given
//   by multiple indices. The default is sum, with multiplication, replacement
//   and omission also being supported. The use of replacement and omission is
//   not recommended, as they are not commutative and depend on the order of the
//   input.
//
// Size: An array with one element per index. Each argument should be a positive
//   integer or NULL. A NULL element is only allowed if the corresponding index
//   is categorical, in which case the cardinality is used instead. This is used
//   to specify the dimensions of the result. The default value is an array full
//   of NULL values.
//
// Check: A boolean specifying whether inputs should be checked to ensure that
//   they fall in the range given by the size array. Not performing a check is
//   faster but riskier, as illegal input will result in a segmentation fault or
//   overwriting random memory in the system. The default value is true.
//
// Error: A boolean specifying whether an error should be thrown when illegal
//   input is given. This argument is only used if check is true. If false,
//   illegal input is simply ignored. The default value is false.
function Build_Fixed_Array(array $t_args, array $inputs, array $outputs)
{
    // Class name is randomly generated.
    $className = generate_name("BuildFixedArray");

    // The number of indices provided.
    $numIndices = count($inputs) - 1;

    // At least two arguments must be provided, as using only one argument would
    // be no better than a simple aggregate on the single input.
    grokit_assert($numIndices >= 1,
                  'Build Array: at least 2 inputs expected.');

    // Processing of inputs is done here.
    $indices = array_slice($inputs, 0, $numIndices, true);

    // Each index is expected to be a scalar. Currently, no scalar property is
    // implemented so each input is simply ensured to be categorical or numeric.
    foreach ($indices as $name => $type)
        grokit_assert($type->is('numeric') || $type->is('categorical'),
                      "Build Array: $name has illegal type $type.");
    $indices = array_values($indices);

    // The dimensions of the data is computed.
    $size = get_default($t_args, 'size', array_fill(0, $numIndices, NULL));
    grokit_assert(is_array($size) && count($size) == $numIndices,
                  'Build Array: size should have one element per index.');
    $size = array_values($size);
    foreach ($size as $key => &$value)
        if (is_null($value))
            if ($indices[$key]->is('categorical'))
                $value = $indices[$key]->get('cardinality');
            else
                grokit_error('Build Array: illegal null for size entry $key.');
        else if (is_integer($value))
            grokit_assert($value > 0, 'Build Array: non-positive size given.');
        else
            grokit_error('Build Array: illegal size type: ' . gettype($value));

    $check = get_default($t_args, 'check', true);
    grokit_assert(is_bool($check), 'Build Array: check should be a boolean.');

    if ($check) {
        $error = get_default($t_args, 'error', false);
        grokit_assert(is_bool($check), 'Build Array: error is not a boolean.');
    } else if (array_key_exists('error', $t_args)) {
        grokit_warning('Build Array: error given when check is false.');
    }

    // The type of the data being accumulated.
    $type = array_get_index($inputs, $numIndices);

    // The structure of the result is computed.
    $size = array_merge($size, get_dimension($type));
    $dim = count($size);
    grokit_assert($dim <= 3,
                  "Build Fixed Array: cannot build array of dimension $dim");
    $types = [NULL, 'Col', 'Mat', 'Cube'];
    $result = $types[$dim]
            . '::fixed<' . collapse(', ', $size) . '>';

    $isScalar = $dim == $numIndices;
?>

using namespace arma;

class <?=$className?>;

class <?=$className?> {
 public:
  // The type of result.
  using Array = <?=$result?>;
}

<?

}
?>
