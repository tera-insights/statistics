// http://www.stat.cmu.edu/~ryantibs/papers/median.pdf
<?php
function Median_Binning_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className = $t_args['className'];
    $dimension = $t_args['dimension'];
    $numBins = $t_args['numBins'];

    $sys_headers  = ['armadillo'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
?>
using namespace arma;
using namespace std;

class <?=$className?>ConstantState{
 public:
  static const constexpr int kDimension = <?=$dimension?>;

 private:
  // The current iteration.
  // TODO: Perhaps implement a max iteration which halts and produces bins for
  // the medians.
  long iteration;

  // Whether to save values for this iteration. Only done during the last
  // iteration.
  bool final_iteration;

  // Vectors used during state 2. They are the lower bound of the upper bin
  // and the upper bound of the lower bin. They are used to compute the least
  // and greatest elements in the respective bins.
  vec::fixed<kDimension> min_2;
  vec::fixed<kDimension> max_2;

  // Vectors used to compute the binning intervals
  vec::fixed<kDimension> min;
  vec::fixed<kDimension> max;
  vec::fixed<kDimension> range;

  // The total number of entries in the data. It is needed to know the exact
  // index needed for the median.
  int count;

  // This is the number of elements to the left of the current binning section
  // for each dimension. This is N_L in the Tibshirani paper. This is needed to
  // compute which element in the sorted array is the median.
  uvec::fixed<kDimension> left_counts;

  // What exact state each attribute is in after pre-computation.
  // 1 : recursive binning is still needed to narrow down the range
  // 2 : the median follows exactly between two bins.
  // 3 : the number of items in the binning range is below the sort threshold
  // 4 : the median has been found for this attribute
  // The algorithm concludes once all attributes are in state 4.
  uvec::fixed<kDimension> state;

  // Used to store the value of the medians once they are found
  vec::fixed<kDimension> medians;

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState()
      : min(),
        max(),
        range(),
        count(0),
        iteration(0),
        left_counts(),
        final_iteration(false),
        state(fill::zeros) {
  }
};
<?php
    return [
        'kind'           => 'RESOURCE',
        'name'           => $className . 'ConstantState',
        'system_headers' => $sys_headers,
        'user_headers'   => $user_headers,
        'lib_headers'    => $lib_headers,
        'libraries'      => $libraries,
    ];
}

function Median_Binning(array $t_args, array $inputs, array $outputs)
{
    // Class name randomly generated
    $className = generate_name("KMC");

    // Initialization of local variables from template arguments
    $numBins = $t_args['number.bins'];
    $sortThreshold = $t_args['sort.threshold'];
    $epsilon = get_default($t_args, 'epsilon', 0.0001);

    // The dimension of the data, i.e. how many elements are in each item.
    $dimension = count($inputs);

    // Indices shift
    $indicesShift = $dimension > 1
        ? implode(', ', range(1, ($numBins + 2) * $dimension, $numBins + 2))
        : "1";

    // Array for generating inline C++ code compactly
    $codeArray = array_combine(array_keys($inputs), range(0, $dimension - 1));

    // Setting output type
    array_set_index($outputs, 0, lookupType("base::JSON"));

    $sys_headers  = ['armadillo', 'limits'];
    $user_headers = [];
    $lib_headers  = [];
    $libraries    = ['armadillo'];
?>

using namespace arma;
using namespace std;

class <?=$className?>;

<?  $constantState = lookupResource("statistics::Median_Binning_Constant_State",
                                    ['className' => $className,
                                     'dimension' => $dimension,
                                     'numBins'   => $numBins]
    ); ?>

class <?=$className?> {
 public:
  static const constexpr int kDimension = <?=$dimension?>;

 private:
  // The typical constant state for an iterable GLA.
  const <?=$constantState?>& constant_state;

  // A vector whose elements will be individually set to contain the current
  // item's data. The data is packaged as a vector due to the reliance on
  // armadillo. This exists as a member variable as opposed to a local variable
  // in order to avoid repeated memory allocation.
  vec::fixed<kDimension> item;

  // These are used for the first iteration to compute the starting statistics.

  // This is the sum of the items, used to compute the means.
  vec::fixed<kDimension> sum;

  // This is the sum of the items squared, used to compute the variances.
  vec::fixed<kDimension> sum_squared;

  // This contains 3 column - max, min, range - with an entry for each input.
  mat::fixed<kDimension, 3> extremities;

  // The number of elements processed, used to compute the means.
  unsigned long count;

  // These are used during the binning process.

  // This is used to store the counts after the first iteration.
  umat::fixed<<?= $numBins + 2 ?>, kDimension> counts;

  // This stores the indices of counts to increment, using counts.elem(indices).
  uvec::fixed<kDimension> indices;

  // This is hard coded by co-generation at compile time. It shifts the vector
  // of indices in order to account for each index being from a different row.
  vec::fixed<kDimension> index_shift;

  // This is set equal to the first item in the data set. It is used to check
  // whether every data row is exactly the same. If this is the case, then the
  // medians must be those unique values. This is needed in case of some values
  // appearing more times than the sort threshold, in which case the binning
  // algorithm cannot converge to under the sort threshold.
  vec::fixed<kDimension> unique_item;

  // This records whether the attributes of the data set each have a single
  // unique value. Due to machine precision issues, exact equality is not
  // checked; instead, the difference must be less than a given epsilon.
  uvec::fixed<kDimension> is_unique;

  // This records whether, for each attribute, there has been a value found
  // within the binning range. Upon such a value being found, it is used to
  // check against for uniqueness while computing the above.
  uvec::fixed<kDimension> needs_unique;

  // Whether, for each attribute, the current item is within the binning range.
  uvec::fixed<kDimension> in_range;

  // These are used for state 2. They store the maximal element in the lower
  // bin and the minimal element in the upper bin, respectively.
  vec::fixed<kDimension> lower_max;
  vec::fixed<kDimension> upper_min;

  // These are used for the state 3.

  // The number of elements in each bin for the final iteration.
  uvec::fixed<kDimension> bin_counts;

  // A matrix to store the relevant data entries in memory for sorting. Not all
  // rows are filled completely. It is only used during the final iteration.
  mat::fixed<kDimension, <?=$sortThreshold?>> data;

 public:
  <?=$className?>(const <?=$constantState?> &state)
      : constant_state(state),
        item(),
        count(0),
        extremities(),
        sum(zeros<vec>(kDimension)),
        sum_squared(zeros<vec>(kDimension)),
        counts(zeros<umat>(<?=($numBins + 2)?>, kDimension)),
        indices(),
        index_shift({<?=$indicesShift?>}),
        bin_counts(zeros<uvec>(kDimension)),
        is_unique(ones<uvec>(kDimension)),
        needs_unique(ones<uvec>(kDimension)),
        data()  {
    extremities.col(0).fill(-numeric_limits<double>::infinity());
    extremities.col(1).fill( numeric_limits<double>::infinity());
    lower_max.fill(-numeric_limits<double>::infinity());
    upper_min.fill( numeric_limits<double>::infinity());
  }

  // During the first stage, the necessary summary statistics are computed.
  // The max and min are not strictly necessarily, but are computed in the rare
  // case that they specify a smaller range, i.e. either mean - stddev < min or
  // mean + stddev > max. This adds little to no computation time as the first
  // iteration is most likely IO bound. This takes one iteration.

  // During the second stage, the data is sorted into bins using a simple
  // linear transform and the number of items above and below the bins is kept
  // track of. In addition, the state checks that not all items within the bin
  // range are exactly the same. If they are, then the iterative step cannot
  // narrow down the range. Rather, the medians are just those unique elements.
  // This takes on average log_(# bins) (# items / sort threshold) iterations.

  // It should be noted that in this algorithm, the last bin contains elements
  // with item == max, i.e. the interval is closed on both sides and not open
  // above like it is in the referenced paper.

  // State 2 is for attributes whose median fell inbetween two elements in
  // different bins (meaning the count of the data was even). This is treated
  // separately for both speed (we only need the exteremities of those bins) and
  // because the binning algorithm can not converge if the two bins are the
  // first and last bins (meaning all that bins inbetween were empty).

  // During state 3, the range of the bins has been narrowed enough so that few
  // enough elements remain in the range. These elements are stored in memory
  // and then sorted in the should iterative step. The median is then grabbed.
  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    <?=inputs_to_vector($inputs)?>
    if (constant_state.iteration == 0) {
      sum += item;
      sum_squared += item % item;
<?  foreach ($codeArray as $name => $counter) { ?>
      if (<?=$name?> > extremities(<?=$counter?>, 0))
        extremities(<?=$counter?>, 0) = <?=$name?>;
      if (<?=$name?> < extremities(<?=$counter?>, 1))
        extremities(<?=$counter?>, 1) = <?=$name?>;
<?  } ?>
      count++;

    } else {
      // Only needed this next part if we need to reduce the binning range.
      if (any(constant_state.state == 0)) {
        // Checking for uniquness first.
<?  foreach ($codeArray as $name => $counter) { ?>
        if (constant_state.state(<?=$counter?>) == 0) {
          if (   needs_unique(<?=$counter?>)
              && item(<?=$counter?>) <= constant_state.max(<?=$counter?>)
              && item(<?=$counter?>) >= constant_state.min(<?=$counter?>)) {
            needs_unique(<?=$counter?>) = 0;
            unique_item(<?=$counter?>) = item(<?=$counter?>);
          }

          if (   is_unique(<?=$counter?>)
              && item(<?=$counter?>) <= constant_state.max(<?=$counter?>)
              && item(<?=$counter?>) >= constant_state.min(<?=$counter?>)
              && !needs_unique(<?=$counter?>)
              && abs(unique_item(<?=$counter?>) - item(<?=$counter?>)) > <?=$epsilon?>) {
            is_unique(<?=$counter?>) = 0;
          }
        }

<?  } ?>
        // These two transformations sort the data into bins (including outside
        // the binning range).
        item -= constant_state.min;
        item /= <?=(1 / $numBins)?> * constant_state.range;

        // This normalizes the data. If data fell outside the binning range, it
        // is placed into the first or last bin. If the normalized index was
        // exactly equal to the number of bins, its index is decremented in
        // order to make so that values equal to the max of the binning range
        // are binned.
        item.transform([](double index) {
            return index < 0
              ? -1
              : index == <?=$numBins?>
                ? <?=$numBins - 1?>
                : index > <?=$numBins?>
                  ? <?=$numBins?>
                  : index;
        });

        // This accounts for the matrix columns starting at 0, # bins, ...
        // It also increments the indices as necessitated by the above
        // transforms, i.e. values were possibly -1.
        item += index_shift;
        indices = conv_to<uvec>::from(item);
        counts.elem(indices) += 1;
      }

      // This represents state 2. The above stuff is not actually needed for
      // attributes in state 2, but vectorized instructions make it pointless to
      // restrict the above operations and extra computations would be needed
      // anyway to construct the item and then choose which columns need a bin
      // count incremented.
<?  foreach ($codeArray as $name => $counter) { ?>
      if (   constant_state.state(<?=$counter?>) == 2
          && (<?=$name?> > lower_max(<?=$counter?>))
          && (<?=$name?> < constant_state.max_2(<?=$counter?>)))
        lower_max(<?=$counter?>) = <?=$name?>;
      if (   constant_state.state(<?=$counter?>) == 2
          && (<?=$name?> < upper_min(<?=$counter?>))
          && (<?=$name?> > constant_state.min_2(<?=$counter?>)))
        upper_min(<?=$counter?>) = <?=$name?>;
<?  } ?>

    // State 3. Elements are stored in memory.
<?  foreach ($codeArray as $name => $counter) { ?>
      if (   constant_state.state(<?=$counter?>) == 3
          && <?=$name?> >= constant_state.min(<?=$counter?>)
          && <?=$name?> <=  constant_state.max(<?=$counter?>)) {
        data(<?=$counter?>, bin_counts(<?=$counter?>)) = <?=$name?>;
        bin_counts(<?=$counter?>)++;
      }
<?  } ?>
    }
  }

  // Each stage of AddState is very simple as most maintained statistics are
  // additive. Other stuff such as uniqueness, upper_min, and lower_max are
  // aggregated in an obvious manner.
  void AddState(<?=$className?> &other) {
    if (constant_state.iteration == 0) {
      sum += other.sum;
      sum_squared += other.sum_squared;
      extremities.col(0) =
          max(join_rows(other.extremities.col(0), extremities.col(0)), 1);
      extremities.col(1) =
          min(join_rows(other.extremities.col(1), extremities.col(1)), 1);
    count += other.count;
    } else {
      if (any(constant_state.state == 0))
        counts += other.counts;
      vec difference = abs(unique_item - other.unique_item);
      for (int counter = 0; counter < kDimension; counter ++) {
        if (constant_state.state(counter) == 0) {
          is_unique(counter) = is_unique(counter)
              && other.is_unique(counter)
              && difference(counter) < <?=$epsilon?>;
        } else if (constant_state.state(counter) == 2) {
          lower_max(counter) = max(lower_max(counter), other.lower_max(counter));
          upper_min(counter) = min(upper_min(counter), other.upper_min(counter));
        } else if (constant_state.state(counter) == 3) {
          if (other.bin_counts(counter) > 0) {
            count = bin_counts(counter);
            data(counter, span(count, count + other.bin_counts(counter) - 1)) =
                other.data(counter, span(0, other.bin_counts(counter) - 1));
            bin_counts(counter) += other.bin_counts(counter);
          }
        }
      }
    }
  }

  // During pre-computation, the constate state is initialized to begin binning.

  // During the first state, the counts are considered and the interval(s) for
  // the median are computed and the constant state min/max are updated. In the
  // rare case that the count of the data is even and the median is the average
  // of elements in 2 different bins, we switch that attribute to state 2.
  // Otherwise, only one interval is kept and hence the range of the median is
  // narrowed by a factor equal to the number of bins.

  // In the third state, the in-memory data is sorted and the median is
  // grabbed based on how many elements fell to the left of the in-memory data.
  bool ShouldIterate(<?=$constantState?>& modible_state) {
    if (constant_state.iteration == 0) {
      vec::fixed<kDimension> mean(sum / count);
      vec::fixed<kDimension> var(sum_squared / count - mean % mean);
      var *= 1.0 + 1.0 / (count - 1);
      vec::fixed<kDimension> std_dev(sqrt(var));
      modible_state.min = max(join_rows(mean - std_dev, extremities.col(1)), 1);
      modible_state.max = min(join_rows(mean + std_dev, extremities.col(0)), 1);
      modible_state.range = modible_state.max - modible_state.min;
      modible_state.count = count;
      modible_state.iteration++;
      return true;

    } else {
      // This is used to check against the sort threshold.
      uvec::fixed<kDimension> median_counts;

      for (int counter = 0; counter < kDimension; counter++) {
        if (modible_state.state(counter) == 0) {
          if (is_unique(counter)) {
            modible_state.state(counter) = 4;
            modible_state.medians(counter) = unique_item(counter);
          } else {
            int total = counts(0, counter);
            int r_counter=  1;
            for (; total < (modible_state.count + 1) / 2; r_counter++) {
              total += counts(r_counter, counter);
            }
            // This shift is to ensure that bin b (starting at 0) contains the
            // range (min + b * increment, min + (b + 1) * increment, where the
            // increment is (range / # bins). It is shifted by two to account
            // for binning starting at index 1 and the extra increment of the
            // for loop.
            r_counter -= 2;
            double increment = modible_state.range(counter) / <?=$numBins?>;

            // The count is even and the median falls exactly between two bins.
            // We need to know what two bins are needed as they are not
            // necessarily contiguous due to empty bins between them.
            if (   modible_state.count % 2 == 0
                && total == modible_state.count / 2) {
              modible_state.state(counter) = 2;
              int second_index = r_counter + 2;
              while (counts(second_index, counter) == 0)
                second_index++;
              modible_state.min(counter) += r_counter * increment;
              modible_state.max(counter) = modible_state.min(counter)
                  + (second_index - r_counter) * increment;
              modible_state.max_2(counter) = modible_state.min(counter)
                  + (r_counter + 1) * increment;
              modible_state.min_2(counter) = modible_state.min(counter)
                  + (second_index - 1) * increment;

            // Otherwise, proceed in an obvious fashion, i.e. keep the one bin
            // that contains the median or the two entries needed to compute it.
            } else {
              modible_state.min(counter) += r_counter * increment;
              modible_state.max(counter) = modible_state.min(counter) + increment;
              if (counts(r_counter + 1, counter) < <?=$sortThreshold?>)
                modible_state.state(counter) = 3;
              modible_state.left_counts(counter) =
                accu(counts(span(0, r_counter), counter));
            }
          }
        } else if (modible_state.state(counter) == 2) {
          modible_state.medians(counter) =
              (lower_max(counter) + upper_min(counter)) / 2;
          modible_state.state(counter) = 4;
        } else if (modible_state.state(counter) == 3) {
          count = modible_state.count;
          rowvec sorted = sort(data(counter, span(0, bin_counts(counter) - 1)));
          double left_count = modible_state.left_counts(counter);
          modible_state.medians(counter) = modible_state.count % 2 == 1
            ? sorted((count - 1) / 2 - left_count)
            : (sorted(count / 2 - left_count)
                + sorted(count / 2 - 1 - left_count)) / 2;
          modible_state.state(counter) = 4;
        }
      }
      modible_state.range = modible_state.max - modible_state.min;
      modible_state.final_iteration = all(median_counts < <?=$sortThreshold?>);
      modible_state.iteration++;
      return any(modible_state.state != 4);
    }
  }

  // The use of "count == 0" is a bit convoluted.

  // In the first stage, it is non zero as it kept track of the count of data.

  // In the second stage, it is incremented once to keep track of whether
  // unique_item was assigned. If all attributes are unique or all attributes
  // are in stage 2a, then it is kept at 1 and ShouldIterate returns false after
  // setting final_iteration to true. Otherwise, count is set to 0. Hence, if
  // another iteration of binning is needed the first branch is used. Also, if
  // stage 2 is concluding (meaning final_iteration is true and ShouldIterate
  // returned false), then again the first branch is used.

  // If the third stage was in effect, it was used in AddStage and hence is not
  // zero. Final iteration is true so the latter branch is used.

  // The branch cannot only be based on final_iteration because that variable is
  // true for the second to last iteration as well as the last. This ensures that
  // the latter branch is indeed only used on the final iteration.
  void GetResult(<?=typed_ref_args($outputs)?>) {
    Json::Value result(Json::objectValue);
    result["iteration"] = (Json::Value::Int64) constant_state.iteration;
    if (all(constant_state.state == 4)) {
      for (int counter = 0; counter < kDimension; counter ++) {
        result["medians"].append(constant_state.medians(counter));
      }
    } else {
      for (int counter = 0; counter < kDimension; counter ++) {
        result["max"].append(constant_state.max(counter));
        result["min"].append(constant_state.min(counter));
      }
    }
    <?=array_keys($outputs)[0]?> = result;
  }
};

<?
    return [
        'kind'            => 'GLA',
        'name'            => $className,
        'system_headers'  => $sys_headers,
        'user_headers'    => $user_headers,
        'lib_headers'     => $lib_headers,
        'libraries'       => $libraries,
        'iterable'        => TRUE,
        'input'           => $inputs,
        'output'          => $outputs,
        'result_type'     => 'single',
        'generated_state' => $constantState,
    ];
}
?>
