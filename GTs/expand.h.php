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

    $inputs_ = $inputs;
    $outputs_ = array_combine(array_keys($inputs_), $outputs);

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
  static const constexpr int kLength = <?=$size?>;

 private:
  // Pointers to the containers.
<?  foreach ($inputs_ as $name => $type) { ?>
  const <?=$type?>* <?=$name?>_ptr;
<?  } ?>

  int index;

 public:
  <?=$className?>() { }

  void ProcessTuple(<?=const_typed_ref_args($inputs_)?>) {
<?  foreach ($inputs_ as $name => $type) { ?>
    <?=$name?>_ptr = &<?=$name?>;
<?  } ?>
    index = 0;
  }

  bool GetNextResult(<?=typed_ref_args($outputs_)?>) {
    if (index == kLength)
      return false;
<?  foreach ($inputs_ as $name => $type) { ?>
    <?=$name?> = (*<?=$name?>_ptr)[index];
<?  } ?>
    index++;
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
