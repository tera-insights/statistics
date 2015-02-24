<?
function Gather(array $t_args, array $inputs, array $outputs) {
    // Class name randomly generated
    $className = generate_name("PreProcess");

    $normal = $t_args['normal'];
    $center = $t_args['center'];
    $scale  = $t_args['scale'];
    $range  = $t_args['range'];

    $inputs_ = array_combine(['vector'], $inputs);
    $count = $inputs_['vector']->get('size');

    $sys_headers  = ['vector'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
    $extra        = [];
    $result_type  = ['single']
?>

using namespace std;
using namespace arma;

class <?=$className?>;

class <?=$className?> {
 public:
  // The length of the input vector.
  static const constexpr unsigned int kCount = <?=$count?>;

 private:
  // Used to compute the range of each attribute.
  vec::fixed<kCount> min;
  vec::fixed<kCount> max;

  // Used to compute the distribution of the data.
  vec::fixed<kCount> mean;
  vec::fixed<kCount> sdev;
  int count;

 public:
  <?=$className?>()
      : count(0) {
    min.fill( datum::inf);
    max.fill(-datum::inf);
    mean.fill(0);
    sdev.fill(0);
  }

  void AddItem(<?=const_typed_ref_args($inputs_)?>) {
<?  if ($center || $scale) { ?>
    mean += vector;
    count++;
<?  } ?>
<?  if ($scale) ?>
    sdev += vector % vector;
<?  if ($scale) { ?>
    min = min(min, vector);
    max = max(max, vector);
<?  } ?>
  }

  void AddState(<?=$className?>& other) {
<?  if ($center || $scale) { ?>
    mean  += other.vector;
    count += other.count;
<?  } ?>
<?  if ($scale) ?>
    sdev  += other.sdev;
<?  if ($range) { ?>
    min = min(min, other.min);
    max = max(max, other.max);
<?  } ?>
  }

  void FinalizeState() {
<?  if ($center || $scale) { ?>
    mean /= count;
<?  } ?>
<?  if ($scale) ?>
    sdev = sqrt(sdvec / count - mean % mean);
<?  if ($range) { ?>
    min = min(min, other.min);
    max = max(max, other.max);
<?  } ?>
  }
};

<?
    return [
        'kind'           => 'GLA',
        'name'           => $className,
        'system_headers' => $sys_headers,
        'user_headers'   => $user_headers,
        'lib_headers'    => $lib_headers,
        'libraries'      => $libraries,
        'extra'          => $extra,
        'iterable'       => false,
        'input'          => $inputs,
        'output'         => $outputs,
        'result_type'    => $result_type,
    ];
}
?>
