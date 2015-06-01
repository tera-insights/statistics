<?
function Expand(array $t_args, array $inputs, array $outputs, array $states)
{
    // Class name is randomly generated.
    $className = generate_name("Expand");

    // Validating inputs and outputs.
    $count = count($inputs);
    grokit_assert($count > 0, 'Expand: received 0 inputs.');
    grokit_assert(in_array(count($outputs), [$count, $count + 1]),
                  'Expand: an invalid number of outputs was received.');

    foreach (array_values($inputs) as $index => $type) {
        grokit_assert($type->is('vector') || $type->is('array'),
                      "Expanded: $type is invalid.");
        if ($index == 0)
            $size = $type->get('size');
        else
            grokit_assert($type->get('size') == $size,
                          'Expand: all inputs must have the same size');
    }

    // Setting output types.
    for ($index = 0; $index < $count; $index++)
        array_set_index($outputs, $index,
                        array_get_index($inputs, $index)->get('type'));

    // Return values.
    $sys_headers  = ['armadillo'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
?>

using namespace arma;

class <?=$className?>;

class <?=$className?> {
 public:
  // The size of each container.
  static const constexpr int kLength = <?=$size?>

 private:
  // The vector to be placed on top of the input array.

 public:
  <?=$className?>(const <?=$constantState?>& state)
      : constant_state(state)  {
  }

  bool ProcessTuple(<?=process_tuple_args($inputs, $outputs)?>) {
<?  if ($normal) { ?>
    item = normalise(rowvec(<?=$point?>.data(), kLength));
    distances = item * constant_state.neighbors;
    distances.max(prediction);
<?  } else { ?>
    item = colvec(<?=$point?>.data(), kLength);
    neighbors.each_col() -= item;
    distances = sum(square(neighbors));
    distances.min(prediction);
<?  } ?>
    <?=$class?> = prediction;
    return true;
  }
};

<?
    return [
        'kind'            => 'GT',
        'name'            => $className,
        'system_headers'  => $sys_headers,
        'user_headers'    => $user_headers,
        'lib_headers'     => $lib_headers,
        'libraries'       => $libraries,
        'input'           => $inputs,
        'output'          => $outputs,
        'result_type'     => 'multi',
    ];
}
?>
