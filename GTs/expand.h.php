<?
function Expand(array $t_args, array $inputs, array $outputs, array $states)
{
    // Class name is randomly generated.
    $className = generate_name("Expand");

    // Processing of inputs.
    $num_inputs = count($inputs);
    grokit_assert($num_inputs > 0, 'Expand: received 0 inputs.');

    foreach (array_values($inputs) as $index => $type) {
        grokit_assert($type->is('vector') || $type->is('array'),
                      "Expanded: $type is invalid.");
        if ($index == 0)
            $size = $type->get('size');
        else
            grokit_assert($type->get('size') == $size,
                          'Expand: all inputs must have the same size');
    }

    $inputs_ = $inputs;

    // Processing of outputs.
    $num_outputs = count($outputs);
    $produce_index = $num_outputs == $num_inputs + 1;
    grokit_assert($produce_index || $num_outputs == $num_inputs,
                  'Expand: an invalid number of outputs was received.');

    for ($index = 0; $index < $num_outputs; $index++)
        array_set_index($outputs, $index,
                        array_get_index($inputs, $index)->get('type'));

    $output_keys = array_keys($inputs_);
    if ($produce_index) {
        array_set_index($outputs, $num_outputs, lookupType('int'));
        $output_keys[] = 'index';
    }
    $outputs_ = array_combine($output_keys, $outputs);

    // Return values.
    $sys_headers  = [];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = [];
?>

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

  // The iterator used to traverse the containers.
  int iter;

 public:
  <?=$className?>() {}

  void ProcessTuple(<?=const_typed_ref_args($inputs_)?>) {
<?  foreach ($inputs_ as $name => $type) { ?>
    <?=$name?>_ptr = &<?=$name?>;
<?  } ?>
    iter = 0;
  }

  bool GetNextResult(<?=typed_ref_args($outputs_)?>) {
    if (iter == kLength)
      return false;
<?  foreach ($inputs_ as $name => $type) { ?>
    <?=$name?> = (*<?=$name?>_ptr)[iter];
<?  } ?>
<?  if ($produce_index) { ?>
    index = iter;
<?  } ?>
    ++iter;
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
