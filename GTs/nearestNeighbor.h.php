<?
function Nearest_Neighbor_Constant_State(array $t_args)
{
    // Initialization of local variables from template arguments.
    $className = $t_args['className'];
    $states    = $t_args['states'];
    $normal    = $t_args['normal'];

    // Values to be used in C++ code.
    $state = array_keys($states)[0];
    $class = $states[$state];

    // Return values.
    $sys_headers = ['armadillo'];
    $user_headers = [];
    $lib_headers = [];
    $libraries = ['armadillo'];
?>

using namespace arma;

class <?=$className?>ConstantState {
 private:
  // The length of each point.
  static const constexpr unsigned int kLength = <?=$class?>::kLength;

  // The number of neighbors.
  static const constexpr unsigned int kCardinality = <?=$class?>::kCardinality;

 public:
  // The matrix containing the neighboring points by column.
  mat::fixed<kLength, kCardinality> neighbors;

  friend class <?=$className?>;

  <?=$className?>ConstantState(<?=const_typed_ref_args($states)?>) {
<?  if ($normal) { ?>
    neighbors = normalise(<?=$state?>.GetCharts());
<?  } else { ?>
    neighbors = <?=$state?>.GetCharts();
<?  } ?>
  }
};
<?
    return [
        'kind'           => 'RESOURCE',
        'name'           => $className . 'ConstantState',
        'system_headers' => $sys_headers,
        'user_headers'   => $user_headers,
        'lib_headers'    => $lib_headers,
        'libraries'      => $libraries,
        'configurable'   => false,
    ];
}

function Nearest_Neighbor(array $t_args, array $inputs, array $outputs, array $states)
{
    // Setting output type
    $state = array_get_index($states, 0);
    array_set_index($outputs, 0, $state->get('key')->lookup());

    // Class name is randomly generated.
    $className = generate_name("NN");

    // Initialization of local variables from template arguments
    $normal = get_default($t_args, 'normalization', false);

    // Initializiation of argument names.
    $point = array_keys($inputs)[0];   // Name of the input point as an array.
    $class  = array_keys($outputs)[0];  // Name of the predicted type outputted.

    // Return values.
    $sys_headers = ['armadillo'];
    $user_headers = [];
    $lib_headers = [];
    $libraries = ['armadillo'];
?>

using namespace arma;

class <?=$className?>;

<?  $constantState = lookupResource(
        "statistics::Nearest_Neighbor_Constant_State",
        ['className' => $className, 'states' => $states, 'normal' => $normal]
    ); ?>

class <?=$className?> {
 public:
  // The length of each point.
  static const constexpr unsigned int kLength = <?=$state?>::kLength;

  // The number of neighbors.
  static const constexpr unsigned int kCardinality = <?=$state?>::kCardinality;

 private:
  // The constant state used to hold matrix..
  const <?=$constantState?>& constant_state;

  // The vector to be placed on top of the input array.
<?  if ($normal) { ?>
  rowvec::fixed<kLength> item;
<?  } else { ?>
  colvec::fixed<kLength> item;
<?  } ?>

  // The matrix used to copy the constant state neighbors.
  mat::fixed<kLength, kCardinality> neighbors;

  // The vector containing the distances to each neighbor.
  rowvec::fixed<kCardinality> distances;

  // The predicted class as an unsigned integer.
  uword prediction;

 public:
  <?=$className?>(const <?=$constantState?>& state)
      : constant_state(state),
        item(),
        neighbors(constant_state.neighbors) {
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
        'result_type'     => 'single',
        'generated_state' => $constantState,
    ];
}
?>
