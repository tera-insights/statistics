 <?php
require_once "grokit_base.php";

function FunctionName(array $t_args, array $inputs, array $outputs)
{
    // Class name randomly generated
    $className = generate_name("ABBR");

?>

using namespace arma;
using namespace std;

class <?=$className?>;

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
        'iterable'         => FALSE,
        'input'            => $inputs,
        'output'           => $outputs,
        'result_type'      => 'single',
    );
}
?>