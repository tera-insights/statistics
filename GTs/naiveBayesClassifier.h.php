 <?php
require_once "grokit_base.php";
// This should only be used for data with at least one numeric attribute. It
// will throw an error otherwise.
function Naive_Bayes_Classifier_Predict_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className = $t_args['className'];
    $numN = $t_args['numN'];
?>
using namespace arma;
using namespace std;

class <?=$className?>ConstantState{
 private:
  // The current iteration which is to be compared to $maxIteration.
  long iteration;

  // Vectors used to compute the binning intervals
  vec::fixed<<?=$numN?>> min;
  vec::fixed<<?=$numN?>> max;
  vec::fixed<<?=$numN?>> range;

  // The total number of entries in the data. It is needed to know the exact
  // index needed for the median.
  int count;

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState()
      : min(),
        max(),
        range(),
        count(0),
        iteration(0) {
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

function Naive_Bayes_Classifier_Predict(array $t_args, array $inputs, array $outputs)
{
    // Setting output type
    array_set_index($outputs, 0, lookupType("base::JSON"));

    // Class name randomly generated
    $className = generate_name("NBC");

    // Initialization of local variables from template arguments

    // This specifies how many standard deviations wide the histogram is in each
    // direction.
    $histogramWidthFactor = $t_args['histogram.width.factor'];

    // The number of bins in the histogram range.
    $numBins = $t_args['number.bins'];

    $numericDistributions = array(
        'gaussian',
        'inverseGaussian',
        'gamma',
        'poisson'
    );

    // Processing the response variable. It must be a categorical input.
    $response = strval($t_args['response']);

    grokit_assert(   array_key_exists($response, $inputs)
                  && $inputs[$response]->is("categorical"),
                  "Response does not specify a categorical input.");

    $numC = $inputs[$response]->get("cardinality");

    // Splitting the inputs into numeric predictors, categorical predictors, and
    // the response.
    $factors = array();
    $numeric = array();

    foreach ($inputs as $name => $type) {
        if ($type->is("numeric"))
            $numeric[] = $name;
        else if ($type->is("categorical") && $name != $response)
            $factors[$name] = $type->get("cardinality");
        else
            grokit_assert($name == $response, "Unusual type encountered.");
    }

    $numF = count($factors);
    if ($hasFactors = ($numF > 0))
           $maxPower = max($factors);
    $numN = count($numeric);

    grokit_assert($numN > 0, "Number of numeric attributes is 0.");

    // Indices shift
    $indicesShiftN = $numN > 1
        ? implode(', ', range(1, ($numBins + 2) * $numN, $numBins + 2))
        : "1";
    $indicesShiftF = $numF > 1
        ? implode(', ', range(0, $maxPower * $numF - 1, $maxPower))
        : "0";
?>

using namespace arma;
using namespace std;

class <?=$className?>;

<?  $constantState = lookupResource(
        'statistics::Naive_Bayes_Classifier_Constant_State',
        array(
            'className' => $className,
            'numN'      => $numN,
        )
    ); ?>

class <?=$className?> {
 private:
  // The typical constant state for an iterable GLA.
  const <?=$constantState?> & constant_state;

  // A vector whose elements will be individually set to contain the current
  // item's data. The data is packaged as a vector due to the reliance on
  // armadillo. This exists as a member variable as opposed to a local variable
  // in order to avoid repeated memory allocation.
  vec::fixed<<?=$numN?>> item;

<?  if ($hasFactors) { ?>
  // We need two different vectors here because numeric and factor attributes
  // are treated quite differently;
  vec::fixed<<?=$numF?>> factors;

<?  } ?>
  // These are used for the first iteration to compute the starting statistics.

  // This is the sum of the items, used to compute the histogram intervals.
  vec::fixed<<?=$numN?>> sum;

  // This is the sum of the items squared, used to compute the variances.
  vec::fixed<<?=$numN?>> sum_squared;

  // This contains 3 column - max, min, range - with an entry for each input.
  mat::fixed<<?=$numN?>, 3> extremities;

  // The number of elements processed, used to compute the means.
  unsigned long count;

  // These are used for the second iteration to count the data.

<?  if ($hasFactors) { ?>
  // This stores the counts of factors. Each slice corresponds to one class of
  // the response and each column to a factor.
  ucube::fixed<<?=$maxPower?>, <?=$numF?>, <?=$numC?>> factor_counts;

  // This is hard coded by co-generation at compile time. It shifts the vector
  // of indices in order to account for each index being from a different row.
  vec::fixed<<?=$numF?>> index_shift_f;

<?  } ?>
  // This stores the histograms of numeric data. Each slice corresponds to one
  // class of the response and each column to the histogram of an attribute.
  ucube::fixed<<?= $numBins + 2 ?>, <?=$numN?>, <?=$numC?>> numeric_counts;

  // This stores the indices of counts to increment, using counts.elem(indices).
  uvec::fixed<<?=$numN?>> indices;

  // This is hard coded by co-generation at compile time. It shifts the vector
  // of indices in order to account for each index being from a different row.
  vec::fixed<<?=$numN?>> index_shift_n;

  // This stores the counts of each class of the response.
  uvec::fixed<<?=$numC?>> class_counts;

 public:
  <?=$className?>(const <?=$constantState?> &state)
      : constant_state(state),
        item(),
        count(0),
        extremities(),
        sum(zeros<vec>(<?=$numN?>)),
        sum_squared(zeros<vec>(<?=$numN?>)),
<?  if ($hasFactors) { ?>
        factors_counts(zeros<uvec>(<?=$maxPower?>, <?=$numF?>, <?=$numC?>)),
        index_shift_f({<?=$indicesShiftF?>}),
<?  } ?>
        numeric_counts(zeros<ucube>(<?=$numBins + 2?>, <?=$numN?>, <?=$numC?>)),
        indices(),
        index_shift_n({<?=$indicesShiftN?>}),
        class_counts(zeros<uvec>(<?=$numC?>)) {
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    <?=inputs_to_vector(array_flip($numeric))?>
    if (constant_state.iteration == 0) {
      sum += item;
      sum_squared += item % item;
<?  foreach ($numeric as $counter => $name) { ?>
      if (<?=$name?> > extremities(<?=$counter?>, 0))
        extremities(<?=$counter?>, 0) = <?=$name?>;
      if (<?=$name?> < extremities(<?=$counter?>, 1))
        extremities(<?=$counter?>, 1) = <?=$name?>;
<?  } ?>
      count++;
    } else {
      // These two transformations sort the data into bins (including outside
      // the binning range).
      item -= constant_state.min;
      item /= <?=(1 / $numBins)?> * constant_state.range;

      // This normalizes the data. If data fell outside the binning range, it is
      // placed into the first or last bin. If the normalized index was exactly
      // equal to the number of bins, its index is decremented in order to make
      // so that values equal to the max of the binning range are binned.
      item.transform([](double index) {
          return index < 0
            ? -1
            : index == <?=$numBins?>
              ? <?=$numBins - 1?>
              : index > <?=$numBins?>
                ? <?=$numBins?>
                : index;
      });

      // This accounts for the matrix columns starting at 0, # bins, ... It also
      // increments the indices as necessitated by the above transforms, i.e.
      // values were possibly -1.
      item += index_shift_n;
      indices = conv_to<uvec>::from(item);
      numeric_counts.slice(<?=$response?>).elem(indices) += 1;

<?  if ($hasFactors) { ?>
      // Very little works needs to be done for categorical data. We are just
      // tabulating each occurrence.
      <?=inputs_to_vector($factors, 'factors')?>
      factors += index_shift_f;
      factor_counts.slice(<?=$response?>).elem(factors) += 1;

<?  } ?>
      class_counts(<?=$response?>)++;
    }
  }

  void AddState(<?=$className?> &other) {
    if (constant_state.iteration == 0) {
      sum += other.sum;
      sum_squared += other.sum_squared;
      extremities.col(0) = max(
          join_rows(other.extremities.col(0), extremities.col(0)), 1);
      extremities.col(1) = min(
          join_rows(other.extremities.col(1), extremities.col(1)), 1);
    count += other.count;
    } else {
      numeric_counts += other.numeric_counts;
<?  if ($hasFactors) { ?>
      factor_counts += other.factor_counts;
<?  } ?>
      class_counts += other.class_counts;
    }
  }

  bool ShouldIterate(<?=$constantState?> &modible_state) {
    if (constant_state.iteration == 0) {
      vec::fixed<<?=$numN?>> mean(sum / count);
      vec::fixed<<?=$numN?>> var(sum_squared / count - mean % mean);
      var *= 1.0 + 1.0 / (count - 1);
      vec::fixed<<?=$numN?>> std_dev(sqrt(var));
      std_dev *= <?=$histogramWidthFactor?>;
      modible_state.min = max(join_rows(mean - std_dev, extremities.col(1)), 1);
      modible_state.max = min(join_rows(mean + std_dev, extremities.col(0)), 1);
      modible_state.range = modible_state.max - modible_state.min;
      modible_state.count = count;
      modible_state.iteration++;
      return true;
    } else {
      // Transformations that would normally appear are below in GetResult to
      // avoid having to have class variables whose only purpose is to store the
      // result from here so that GetResult can use it.
      modible_state.iteration++;
      return false;
    }
  }

  void GetResult(<?=typed_ref_args($outputs)?>) {
    if (constant_state.iteration > 1) {
      cube numeric_probabilities(<?= $numBins + 2 ?>, <?=$numN?>, <?=$numC?>);
      ucube::iterator first_iter = numeric_counts.begin();
      cube::iterator second_iter = numeric_probabilities.begin();
      ucube::iterator end = numeric_counts.end();
      for(; first_iter != end; ++first_iter) {
        *second_iter = (double) *first_iter;
        ++second_iter;
      }
      numeric_probabilities /= constant_state.count;
<?  if ($hasFactors) { ?>
      cube factor_probabilities(<?=$maxPower?>, <?=$numF?>, <?=$numC?>);
      first_iter = factor_counts.begin();
      second_iter = factor_probabilities.begin();
      end = factor_counts.end();
      for(first_iter; first_iter != end; ++first_iter) {
        *second_iter = (double) *first_iter;
        ++second_iter;
      }
      factor_probabilities /= constant_state.count;
<?  } ?>
      vec class_probabilities(<?=$numC?>);
      for (int counter = 0; counter < <?=$numC?>; counter++)
        class_probabilities(counter) = (double) class_counts(counter);
      class_probabilities /= constant_state.count;
      Json::Value result(Json::objectValue);
      for (int counter = 0; counter < <?=$numN?>; counter ++) {
        result["max"].append(constant_state.max(counter));
        result["min"].append(constant_state.min(counter));
      }
<?  foreach ($numeric as $counter => $name) { ?>
      for (int class_counter = 0; class_counter < <?=$numC?>; class_counter++) {
        for (int counter = 0; counter < <?= $numBins + 2 ?>; counter++) {
          result["numeric probabilities"][class_counter]["<?=$name?>"].append(
              numeric_probabilities(counter, <?=$counter?>, class_counter)
          );
        }
      }
<?  } ?>
<?  if ($hasFactors) { ?>
<?      foreach ($numeric as $counter => $name) { ?>
      for (int class_counter = 0; class_counter < <?=$numC?>; class_counter++) {
        for (int counter = 0; counter < <?= $numBins + 2 ?>; counter++) {
          result["categorical"][class_counter]["<?=$name?>"].append(
              factor_probabilities(counter, <?=$counter?>, class_counter);
          )
        }
      }
<?      } ?>
<?  } ?>
      for (int counter = 0; counter < <?=$numC?>; counter++)
        result["responseProbability"].append(class_probabilities(counter));
      result["count"] = constant_state.count;
<?  $describer = $inputs[$response]->describer('json');  ?>
<?  $describer('result["responseInfo"]'); ?>
<?= ProduceResult(array_keys($outputs)[0], "result") ?>
    }
  }
};

<?
    return array(
        'kind'             => 'GT',
        'name'             => $className,
        'system_headers'   => array('armadillo', 'string', 'iostream',),
        'user_headers'     => array(),
        'iterable'         => TRUE,
        'input'            => $inputs,
        'output'           => $outputs,
        'result_type'      => 'single',
        'generated_states' => array($constantState),
    );
}
?>
