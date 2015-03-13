<?
//  Copyright 2014 Tera Insights, LLC. All Rights Reserved.

// Converts a given number of inputs into a fixed size vector.

// MakeVector uses a hierarchy of types to determine the overall type of the
// vector. The following lists common types in increasing order.
// 1. 8-bit integer value.
// 2. 16-bit integer value.
// 3. Single precision floating point.
// 4. 32-bit integer value.
// 5. Double precision floating point.
// 6. 64-bit integer value.
// 7. Double extended precision floating point.
// If both floating point and integral types are present, the next highest
// floating point type is used. For example, if both (3) and (4) are seen,
// then (5), at least, would be used as the overall type.

// Currently (7) is not implemented. If both a real type and a 64-bit integer
// are used, an error is thrown.

function MakeVector(array $inputs, array $t_args) {
    $count = count($inputs);
    $direction = get_default($t_args, 'direction', 'col');
    $type = get_default($t_args, 'type', lookupType('base::double'));

    grokit_warning_assert($count, 'MakeVector: 0 inputs received.');

    // Processing of input types.
    $size = 0;
    $maxRealSize = $maxIntegralSize = 0;
    $typeErrorMessage = '';
    foreach ($inputs as $counter => $input) {
        if ($input->is('vector') || $input->is('array'))
            $size += $input->get('size');
        else if ($inputs->is('categorical') || $inputs->is('numeric'))
            $size++;
        else
            $typeErrorMessage .= "Input [$counter] has type $input.";
        /* if ($input->is('real')) */
        /*     $maxRealSize = max($maxRealSize, $input->get('size.bytes')); */
        /* else */
        /*     $maxIntegralSize = max($maxIntegralSize, $input->get('size.bytes')); */
    }
    /* $maxSize = max($maxRealSize, $maxIntegralSize); */

    grokit_assert($typeErrorMessage == '',
                  "MakeVector: Non-numeric inputs:\n" . $typeErrorMessage);

    /* grokit_assert(is_null($type) || $type->get('size.bytes') >= $maxSize, */
    /*               'The given vector type is insufficient to hold the inputs.'); */

    /* if (is_null($type)) { */
    /*     // Null padding used to reduce case-work when choosing type. */
    /*     $typeArray = [ */
    /*         null, lookupType('base::byte'), */
    /*         null, lookupType('base::smallint'), */
    /*         lookupType('base::float'), lookupType('base::int'), */
    /*         lookupType('base::double'), lookupType('base::bigint'), */
    /*     ]; */
    /*     grokit_assert($maxIntegralSize < 8 || $maxRealSize == 0, */
    /*                    'MakeVector: extended double not currently supported.'); */
    /*     if ($maxRealSize == $maxSize && $maxIntegralSize < $maxSize) */
    /*         $index = 2 * log($maxRealSize, 2); */
    /*     else if ($maxRealSize > 0 && $maxIntegralSize == $maxSize) */
    /*         $index = 2 * log($maxIntegralSize, 2) + 2; */
    /*     else */
    /*         $index = 2 * log($maxIntegralSize, 2) + 1; */
    /*     $type = $typeArray[$index]; */
    /* } */

    foreach (range(0, $count - 1) as $counter)
      $arguments[] = 'arg' . $counter;

    $inputs_ = array_combine($arguments, $inputs);

    grokit_assert(is_string($direction) && in_array($direction, ['row', 'col']),
                  'MakeVector: [direction] argument must be "row" or "col".');
    grokit_assert(is_datatype($type) && $type->is('numeric'),
                  'MakeVector: [type] argument must be a numeric datatype.');

    $output = lookupType(
        'statistics::fixed_vector',
        ['direction' => $direction, 'size' => $size, 'type' => $type,
         'inputs' => array_values($inputs)]
    );

    $funName = generate_name('MakeVector');
    $sys_headers  = ['armadillo'];
    $user_headers = [];
    $lib_headers  = [];
?>

<?=$output?> <?=$funName?>(<?=const_typed_ref_args($inputs_)?>) {
  <?=$output?> result;
<?  $count = 0;
    foreach ($inputs_ as $arg => $type) {
        if ($type->is('categorical') || $type->is('numeric')) { ?>
  result(<?=$count?>) = <?=$arg?>;
<?          $count++;
        } else {
            $end = $count + $type->get('size') - 1;
            if ($type->is('array')) { ?>
  result.subvec(<?=$count?>, <?=$end?>) =
      <?=$output?>(<?=$arg?>.data(), <?=$type->get('size')?>);
<?          } else { ?>
  result.subvec(<?=$count?>, <?=$end?>) = <?=$arg?>;
<?          }
            $count = $end + 1;
        }
    } ?>
  return result;
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

