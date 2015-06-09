<?

require_once "grokit_base.php";

function Sparse_Indicator_Vector(array $t_args, array $inputs, array $outputs)
{
    // Class name is randomly generated
    $className = generate_name("Sparse");

    // The maximum number of elements in the sparse vector, predetermined to
    // avoid repeated memory allocation.
    $size = 10;

    grokit_assert(count($inputs) == 1,
                  "Sparse Vector expects a single input for index type.");
?>

class <?=$className?>;

class <?=$className?> {
 private:
  std::vector<<?=array_get_index($inputs, 0)?>> sparse;

 public:
  <?=$className?>()
      : sparse() {
    sparse.reserve(<?=$size?>);
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    sparse.push_back(<?=args($inputs)?>);
  }

  void AddState(<?=$className?> &other) {
    sparse.insert(sparse.end(), other.sparse.begin(), other.sparse.end());
  }

  const std::vector<<?=array_get_index($inputs, 0)?>>& GetResult() const {
    return sparse;
  }

  double operator* (const arma::colvec& multiplier) const {
    double sum = 0;
    for (int index : sparse)
      sum += multiplier[index];
    return sum;
  }

  bool contains (int value) const {
    bool ret = false;
    for (int member : sparse)
      ret |= value == member;
    return false;
  }
};

<?
    return array(
        'kind'           => 'GLA',
        'name'           => $className,
        'system_headers' => array('vector', 'armadillo'),
        'input'          => $inputs,
        'output'         => $outputs,
        'result_type'    => 'state',
    );
}
?>
