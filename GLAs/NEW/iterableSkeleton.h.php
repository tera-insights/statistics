 <?php
require_once "grokit_base.php";

function FunctionName_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className = $t_args['className'];
?>
using namespace arma;
using namespace std;

class <?=$className?>ConstantState{
 private:

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState()
      : {
  }
};
<?php
    return array(
        'kind'           => 'RESOURCE',
        'name'           => $className . 'ConstantState',
        'system_headers' => array('armadillo'),
        'user_headers'   => array(),
    );
}

function FunctionName(array $t_args, array $inputs, array $outputs)
{
    // Class name randomly generated
    $className = generate_name("ABBR");

?>

using namespace arma;
using namespace std;

class <?=$className?>;

<?  $constantState = lookupResource(
        // Function Name Here
        array(
            'className'  => $className,
        )
    ); ?>

class <?=$className?> {
 private:
  // The typical constant state for an iterable GLA.
  const <?=$constantState?> & constant_state;

 public:
  <?=$className?>(const <?=$constantState?> &state)
      : {
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
  }

  void AddState(<?=$className?> &other) {
  }

  bool ShouldIterate(<?=$constantState?> & modible_state) {
  }

  void GetResult(<?=typed_ref_args($outputs)?>) {
};

<?
    return array(
        'kind'             => 'GLA',
        'name'             => $className,
        'system_headers'   => array('armadillo', 'string', 'iostream',
                                    'limits'),
        'user_headers'     => array(),
        'iterable'         => TRUE,
        'input'            => $inputs,
        'output'           => $outputs,
        'result_type'      => 'single',
        'generated_states' => array($constantState),
    );
}
?>