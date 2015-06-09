<?php
require_once "grokit_base.php";

function Median_Binning_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className = $t_args['className'];
    $dimension = $t_args['dimension'];
    $numberBins = $t_args['numberBins'];
?>
using namespace arma;
using namespace std;

class <?=$className?>ConstantState{
 private:
  // The current iteration which is to be compared to $maxIteration.
  long iteration;

  // Whether to save values for this iteration. Only done during the last
  // iteration.
  bool final_iteration;

  // Vectors and scalars  used to compute the binning intervals
  vec::fixed<<?=$dimension?>> min;
  vec::fixed<<?=$dimension?>> max;
  vec::fixed<<?=$dimension?>> range;

  // The total number of entries in the data. It is needed to know the exact
  // index needed for the median.
  int count;

  // This is the number of elements to the left of the current binning section
  // for each dimension. This is N_L in the Tibshirani paper.
  uvec::fixed<<?=$dimension?>> left_counts;

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState()
      : min(),
        max(),
        range(),
        count(0),
        iteration(0),
        left_counts(),
        final_iteration(false) {
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

function Median_Binning(array $t_args, array $inputs, array $outputs)
{
    // Class name is randomly generated
    $className = generate_name("KMC");

    // Initialization of local variables from template arguments

    $numberBins = $t_args['number.bins'];
    $sortThreshold = $t_args['sort.threshold'];

    // The dimension of the data, i.e. how many elements are in each item.
    $dimension = count($inputs);

    // Indices shift
    $indicesShift = implode(', ', range(1, $numberBins * $dimension, $numberBins));

    // Array for generating inline C++ code;
    $codeArray = array_combine(array_keys($inputs), range(0, $dimension - 1));
?>

using namespace arma;
using namespace std;

class <?=$className?>;

<?  $constantState = lookupResource(
        "statistics::Median_Binning_Constant_State",
        array(
            'className'  => $className,
            'dimension'  => $dimension,
            'numberBins' => $numberBins,
        )
    ); ?>

class <?=$className?> {
 private:
  // The typical constant state for an iterable GLA.
  const <?=$constantState?> & constant_state;

  // These are used for the first iteration to compute the starting statistics.

  // A vector whose elements will be individually set to contain the current
  // item's data. The data is packaged as a vector due to the reliance on
  // armadillo. This exists as a member variable as opposed to a local variable
  // in order to avoid repeated memory allocation.
  vec::fixed<<?=$dimension?>> item;

  // This is the sum of the items, used to compute the means.
  vec::fixed<<?=$dimension?>> sum;

  // This is the sum of the items squared, used to compute the variances.
  vec::fixed<<?=$dimension?>> sum_squared;

  // This contains 3 column - max, min, range - with an entry for each input.
  mat::fixed<<?=$dimension?>, 3> extremities;

  // The number of elements processed, used to compute the means.
  long count;

  // These are used during the binning process.

  // This is used to store the counts after the first iteration.
  umat::fixed<<?=$dimension?>, <?=($numberBins + 2)?>> counts;

  // This stores the indices of counts to increment, using counts.elem(indices).
  uvec::fixed<<?=$dimension?>> indices;

  // This stores a vector of 1s. It is created here in order to avoid repeated
  // memory allocation. Calling ones<vec>() repeatedly would be very costly.
  // TODO: Make sure this is actually needed.
  vec::fixed<<?=$dimension?>> increment;

  // This is hard coded by co-generation at compile time. It shifts the vector
  // of indices in order to account for each index being from a different row.
  uvec::fixed<<?=$dimension?>> index_shift;

  // These are used for the final iteration.

  // The number of elements in each bin for the final iteration.
  uvec::fixed<<?=$dimension?>> bin_counts;

  // A matrix to store the relevant data entries in memory for sorting. Not all
  // rows are filled completely. It is only used during the final iteration.
  mat::fixed<<?=$dimension?>, <?=$sortThreshold?>> data;

  // This stores the medians of the data in order to use them in GetResult after
  // computing them in ShouldIterate.
  vec::fixed<<?=$dimension?>> medians;

 public:
  <?=$className?>(const <?=$constantState?> &state)
      : constant_state(state),
        item(),
        count(0),
        extremities(),
        sum(zeros<vec>(<?=$dimension?>)),
        sum_squared(zeros<vec>(<?=$dimension?>)),
        counts(zeros<umat>(<?=$dimension?>, <?=($numberBins + 2)?>)),
        indices(),
        increment(ones<vec>(<?=$dimension?>)),
        index_shift({<?=$indicesShift?>}),
        bin_counts(zeros<uvec>(<?=$dimension?>)),
        data(),
        medians() {
    //extremities.col(0).fill(-numeric_limits<double>::infinity());
    //extremities.col(1).fill( numeric_limits<double>::infinity());
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    <?=inputs_to_vector($inputs)?>
    if (constant_state.iteration == 0) {
<?  foreach ($codeArray as $name => $counter) { ?>
      item[<?=$counter?>] = <?=$name?>;
<?  } ?>
      sum += item;
      sum_squared += item % item;
<?  foreach ($codeArray as $name => $counter) { ?>
      if (<?=$name?> > extremities(<?=$counter?>, 0))
        extremities(<?=$counter?>, 0) = <?=$name?>;
      if (<?=$name?> < extremities(<?=$counter?>, 1))
        extremities(<?=$counter?>, 1) = <?=$name?>;
<?  } ?>
      count++;
    } else if (!constant_state.final_iteration) {
      item -= constant_state.min;
      item /= <?=(1 / $numberBins)?> * constant_state.range;
      indices = conv_to<uvec>::from(item);
      indices.transform([](uword index) {
          return index < 0
              ? -1
              : index > <?=$numberBins?>
                  ? <?=$numberBins?>
                  : index;
      });
      indices += index_shift;
      counts.elem(indices) += ones<uvec>(<?=$dimension?>);
    } else {
<?  foreach ($codeArray as $name => $counter) { ?>
      if (   <?=$name?> >= constant_state.min(<?=$counter?>)
          && <?=$name?> <  constant_state.max(<?=$counter?>)) {
        data(<?=$counter?>, bin_counts(<?=$counter?>)) = <?=$name?>;
        bin_counts(<?=$counter?>)++;
      }
<?  } ?>
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
    } else if (!constant_state.final_iteration) {
      counts += other.counts;
    } else {
      for (int counter = 0; counter < <?=$dimension?>; counter ++) {
        count = bin_counts(counter);
        data(counter, span(count, count + other.count - 1)) =
          other.data(counter, span(0, other.count - 1));
        bin_counts(counter) += other.bin_counts(counter);
      }
    }
  }

  bool ShouldIterate(<?=$constantState?> & modible_state) {
    if (constant_state.iteration == 0) {
      vec::fixed<<?=$dimension?>> mean(sum / count);
      vec::fixed<<?=$dimension?>> var(sum_squared / count - mean % mean);
      var *= 1.0 + 1.0 / (count - 1);
      vec::fixed<<?=$dimension?>> std_dev(sqrt(var));
      extremities.col(2) = extremities.col(0) - extremities.col(1);
      modible_state.min = mean - std_dev;
      modible_state.max = mean + std_dev;
      modible_state.range = 2 * std_dev;
      modible_state.count = count;
      modible_state.iteration++;
      return true;
    } else if (!modible_state.final_iteration) {
      uvec::fixed<<?=$dimension?>> median_counts;
      for (int counter = 0; counter < <?=$dimension?>; counter++) {
        int total = counts(counter, 0);
        int r_counter=  1;
        for (r_counter; total < (modible_state.count + 1) / 2; r_counter++)
          total += counts(counter, r_counter);
        // This shift is to ensure that bin b (starting at 0) contains the range
        // (min + b * range / total bins, min + (b + 1) * range / total_bins).
        r_counter -= 2;
        double increment = modible_state.range(counter) / <?=$numberBins?>;
        // The count is even and the median falls exactly between two bins.
        if (modible_state.count % 2 == 0 && total == modible_state.count / 2) {
          modible_state.min(counter) += r_counter * increment;
          modible_state.max(counter) = modible_state.min(counter)
                                     + 2 * increment;
          cout << counts << endl;
          cout << counts(counter, span(r_counter + 1, r_counter + 2)) << endl;
          cout << sum(counts(counter, span(r_counter + 1, r_counter + 2))) << endl;
          //uword dummy = (uword) sum(counts(counter, span(r_counter + 1, r_counter + 2)), 1);
          //median_counts(counter) = dummy;
              //(uword) sum(counts(counter, span(r_counter + 1, r_counter + 2)));
          //modible_state.left_counts(counter) =
              //(uword) sum(counts(counter, span(0, r_counter)));
        } else {
          modible_state.min(counter) += r_counter * increment;
          modible_state.max(counter) = modible_state.min(counter) + increment;
          cout << "counts" << endl << counts << endl;
          urowvec a = counts(counter, span(r_counter + 1, r_counter + 2));
          cout << "sub" << endl << a << endl;
          cout << "sum" << endl << sum(a) << endl;
          cout << "here" << endl;
          //median_counts(counter) = counts(counter, r_counter + 1);
          //modible_state.left_counts(counter) = (uword) 1;
              //(uword) sum(counts(counter, span(0, r_counter)));
        }
      }
      modible_state.range = modible_state.max - modible_state.min;
      modible_state.final_iteration = all(median_counts < <?=$sortThreshold?>);
      modible_state.iteration++;
      return true;
    } else {
      for (int counter = 0; counter < <?=$dimension?>; counter++) {
        count = bin_counts(counter);
        vec sorted = sort(data(counter, span(0, counter - 1)));
        medians(counter) = modible_state.count % 2 == 1
          ? sorted((count - 1) / 2)
          : (sorted(count / 2) + sorted(count / 2 - 1)) / 2;
      }
      modible_state.iteration++;
      return false;
    }
  }

  void GetResult(<?=typed_ref_args($outputs)?>) {
    Json::Value result(Json::objectValue);
    result["iteration"] = (Json::Value::Int64) constant_state.iteration;
    if (!constant_state.final_iteration || count == 0) {
      for (int counter = 0; counter < <?=$dimension?>; counter ++) {
        result["max"].append(constant_state.max(counter));
        result["min"].append(constant_state.min(counter));
      }
    } else {
      for (int counter = 0; counter < <?=$dimension?>; counter ++) {
        result["medians"].append(medians(counter));
      }
    }
    cout << result << endl;
  }
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
