// http://www.stat.cmu.edu/~ryantibs/papers/median.pdf
 <?php
require_once "grokit_base.php";

function Median_Binning_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className = $t_args['className'];
    $dimension = $t_args['dimension'];
    $numBins = $t_args['numBins'];
?>
using namespace arma;
using namespace std;

class <?=$className?>ConstantState{
 private:
  // The current iteration.
  // TODO: Perhaps implement a max iteration which halts and produces bins for
  // the medians.
  long iteration;

  // Whether to save values for this iteration. Only done during the last
  // iteration.
  bool final_iteration;

  // Whether to modify stage 2 of AddItem. This is needed in the rare case that
  // the count is even and the two elements needed to compute the median are in
  // the first and last bin, respectively. It is a vector in order to account
  // for each attribute separately. Unlike uniqueness, this needs to be recorded
  // on a per attribute basis or a single attribute could hold up the whole
  // process. The values are booleans in the form of unsigned ints.
  uvec::fixed<<?=$dimension?>> stage_2a;

  // Vectors used during stage 2a
  vec::fixed<<?=$dimension?>> min_2a;
  vec::fixed<<?=$dimension?>> max_2a;

  // Vectors used to compute the binning intervals
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
        final_iteration(false),
        stage_2a(zeros<uvec>(<?=$dimension?>)) {
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
    // Class name randomly generated
    $className = generate_name("KMC");

    // Initialization of local variables from template arguments

    $numBins = $t_args['number.bins'];
    $sortThreshold = $t_args['sort.threshold'];
    $quantiles = $t_args['sort.threshold'];

    if (is_int($quantile))
        $quantiles = array($quantiles);

    // TODO: Maybe make non-integer quantiles, refine numBins restriction.

    // This is used to ensure that the quantiles picked will line up with
    // the binning intervals.
    grokit_assert($numBins % 100 == 0,
                  "Number of bins must be a multiple of 100.");

    // Checks to ensure correct format of quantiles.
    grokit_assert(is_array($quantiles), "Quantiles must be in an array.");
    grokit_assert(count($quantiles) > 0, "No quantiles specified.");
    foreach ($quantiles as $quantile)
        grokit_assert(is_int($qauntile), "$quantile is not an integer.");

    // This affects how large the counting datastructure is.
    $numQ = count($quantiles);

    // The dimension of the data, i.e. how many elements are in each item.
    $dimension = count($inputs);

    // Indices shift
    $indicesShift = $dimension > 1
        ? implode(', ', range(1, ($numBins + 2) * $dimension, $numBins + 2))
        : "1";

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
            'numBins' => $numBins,
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
  vec::fixed<<?=$dimension?>> item;

  // These are used for the first iteration to compute the starting statistics.

  // This is the sum of the items, used to compute the means.
  vec::fixed<<?=$dimension?>> sum;

  // This is the sum of the items squared, used to compute the variances.
  vec::fixed<<?=$dimension?>> sum_squared;

  // This contains 3 column - max, min, range - with an entry for each input.
  mat::fixed<<?=$dimension?>, 3> extremities;

  // The number of elements processed, used to compute the means.
  unsigned long count;

  // These are used during the binning process.

  // This data structure is used for different purposes during the second and
  // third iteration.

  // Second iteration:

  // This is used to store the counts during the second iteration. This is used
  // because during the second iteration, we know nothing about how quantiles
  // are distributed and it does not make sense to use an interval map yet. Once
  // we do, we switch over so that each quantile has its own binning range.

  // Third iteration:

  // In each column, there are $numQ separate binning ranges, separated by one
  // interval each, with two bins on the end for a total number of bins equal to
  // $numQ * $numBins + ($numQ - 1) + 2. The binning range for each quantile are
  // guaranteed to be disjoint and ordered due to the second iteration.
  umat::fixed<<?= $numQ * ($numBins + 1) + 1?>, <?=$dimension?>> counts;

  // This stores the indices of counts to increment, using counts.elem(indices).
  uvec::fixed<<?=$dimension?>> indices;

  // This is hard coded by co-generation at compile time. It shifts the vector
  // of indices in order to account for each index being from a different row.
  vec::fixed<<?=$dimension?>> index_shift;

  // This is set equal to the first item in the data set. It is used to check
  // whether every data row is exactly the same. If this is the case, then the
  // medians must be those unique values. This is needed in case of some values
  // appearing more times than the sort threshold, in which case the binning
  // algorithm cannot converge to under the sort threshold.
  vec::fixed<<?=$dimension?>> unique_item;

  // This records whether all attributes of the data set each have a single
  // unique value. Rather than individual uniqueness being considered, this is
  // used to simplify the algorithm. Run time is not greatly affected due to
  // vectorized instructions.
  bool uniqueness;

  // These are used for stage 2a. They store the maximal element in the lower
  // bin and the minimal element in the upper bin, respectively.
  vec::fixed<<?=$dimension?>> lower_max;
  vec::fixed<<?=$dimension?>> upper_min;

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
        counts(zeros<umat>(<?=($numBins + 2)?>, <?=$dimension?>)),
        indices(),
        index_shift({<?=$indicesShift?>}),
        bin_counts(zeros<uvec>(<?=$dimension?>)),
        uniqueness(true),
        data(),
        medians() {
    extremities.col(0).fill(-numeric_limits<double>::infinity());
    extremities.col(1).fill( numeric_limits<double>::infinity());
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
  // Uniqueness is characterized by all attributes having a single respective
  // value, rather than individually. This simplifies the algorithm and doesn't
  // affect run time much due to vectorized instructions. This takes at most
  // log_(# bins) (# items / sort threshold) iterations.

  // It should be noted that in this algorithm, the last bin contains elements
  // with item == max, i.e. the interval is closed on both sides and not open
  // above like it is in the referenced paper.

  // Additionally, in stage 2 there is a substage 2a. It is referred to as a
  // substage as some attributes can be in stage 2a and some in just stage 2.
  // Stage 2a is for attributes whose median fell inbetween two elements in
  // different bins (meaning the count of the data was even). This is treated
  // separately for both speed (we only need the exteremities of those bins) and
  // because the binning algorithm will not converge if the two bins are the
  // first and last bins (meaning all bins inbetween were empty).

  // During the third stage, the range of the bins has been narrowed enough so
  // that few elements remain in the range. These elements are stored in memory
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

    } else if (!constant_state.final_iteration) {
      if (   count == 0
          && all(item >= constant_state.min)
          && all(item <= constant_state.max)) {
        unique_item = item;
        count++;
      }

      if (   uniqueness
          && any(unique_item != item)
          && all(item >= constant_state.min)
          && all(item <= constant_state.max))
        uniqueness = false;

      // Currently performing the second iteration.
      if (constant_state.iteration == 1) {
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
        item += index_shift;
        indices = conv_to<uvec>::from(item);
        counts.elem(indices) += 1;
      } else {
        for (int counter = 0; counter < <?=$dimension?>; counter++) {
          
      }

      // This represents stage 2a. The above stage 2 is not actually needed for
      // attributes in stage 2a, but vectorized instructions make it pointless
      // to restrict the above operations and extra computations would be needed
      // anyway to construct the item and then choose which columns need a bin
      // count incremented.
<?  foreach ($codeArray as $name => $counter) { ?>
      if (   constant_state.stage_2a(<?=$counter?>)
          && (<?=$name?> > lower_max(<?=$counter?>))
          && (<?=$name?> < constant_state.max_2a(<?=$counter?>)))
        lower_max(<?=$counter?>) = <?=$name?>;
      if (   constant_state.stage_2a(<?=$counter?>)
          && (<?=$name?> < upper_min(<?=$counter?>))
          && (<?=$name?> > constant_state.min_2a(<?=$counter?>)))
        upper_min(<?=$counter?>) = <?=$name?>;
<?  } ?>

    // Stage 3. Only attributes who were not on stage 2a are needed.
    } else {
<?  foreach ($codeArray as $name => $counter) { ?>
      if (   !constant_state.stage_2a(<?=$counter?>)
          && <?=$name?> >= constant_state.min(<?=$counter?>)
          && <?=$name?> <=  constant_state.max(<?=$counter?>)) {
        data(<?=$counter?>, bin_counts(<?=$counter?>)) = <?=$name?>;
        bin_counts(<?=$counter?>)++;
      }
<?  } ?>
    }
  }

  // Each stage of AddState is very simple as all maintained statistics are
  // additive. The only other stuff of note is that uniqueness is tracked in an
  // obvious fashion and the in-memory data is aggregated.
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
      uniqueness =    uniqueness && all(unique_item == other.unique_item)
                   && other.uniqueness;
    } else {
      for (int counter = 0; counter < <?=$dimension?>; counter ++) {
        count = bin_counts(counter);
        data(counter, span(count, count + other.count - 1)) =
          other.data(counter, span(0, other.count - 1));
        bin_counts(counter) += other.bin_counts(counter);
      }
    }
  }

  // During the first stage, the constate state is initialized to begin binning.

  // During the second stage, the counts are considered and the interval(s) for
  // the median are computed and the constant state min/max are updated. In the
  // rare case that the count of the data is even and the median is the average
  // of elements in 2 different bins, we switch that attribute to stage 2a.
  // Otherwise, only one interval is kept and hence the range of the median is
  // narrowed by a factor equal to the number of bins.

  // During the third stage, the in-memory data is sorted and the median is
  // grabbed based on how many elements fell to the left of the in-memory data.
  bool ShouldIterate(<?=$constantState?> & modible_state) {
    if (constant_state.iteration == 0) {
      vec::fixed<<?=$dimension?>> mean(sum / count);
      vec::fixed<<?=$dimension?>> var(sum_squared / count - mean % mean);
      var *= 1.0 + 1.0 / (count - 1);
      vec::fixed<<?=$dimension?>> std_dev(sqrt(var));
      modible_state.min = max(join_rows(mean - std_dev, extremities.col(1)), 1);
      modible_state.max = min(join_rows(mean + std_dev, extremities.col(0)), 1);
      modible_state.range = modible_state.max - modible_state.min;
      modible_state.count = count;
      modible_state.iteration++;
      return true;

    } else if (!modible_state.final_iteration) {
      // Only one distinct value for each attribute that might be the median, so
      // those distinct values must be the median.
      if (uniqueness) {
        medians = unique_item;
        modible_state.final_iteration = true;
        modible_state.iteration++;
        return false;
      }

      // All attributes were in stage 2a. Hence, we know what their medians are
      // and we can stop.
      if (all(modible_state.stage_2a)) {
        medians = (modible_state.min_2a + modible_state.max_2a) / 2;
        modible_state.final_iteration = true;
        modible_state.iteration++;
        return false;
      }

      // Otherwise, at least one attributes is neither unique nor in stage 2a,
      // so another iteration of stage 2 is required:

      // This is modified for reasons describe in GetResult.
      count = 0;

      // This is used to check against the sort threshold.
      uvec::fixed<<?=$dimension?>> median_counts;

      for (int counter = 0; counter < <?=$dimension?>; counter++) {
        if (modible_state.stage_2a(counter)) {
          median_counts(counter) = 0;
          modible_state.max_2a(counter) = lower_max(counter);
          modible_state.min_2a(counter) = upper_min(counter);
        } else {
          int total = counts(0, counter);
          int r_counter=  1;
          for (; total < (modible_state.count + 1) / 2; r_counter++) {
            total += counts(r_counter, counter);
          }
          // This shift is to ensure that bin b (starting at 0) contains the range
          // (min + b * range / total bins, min + (b + 1) * range / total_bins).
          // It is shifted by two to account binning starting at 1 and the extra
          // increment of the for loop.
          r_counter -= 2;
          double increment = modible_state.range(counter) / <?=$numBins?>;

          // The count is even and the median falls exactly between two bins. We
          // need to know what two bins are needed as they are not necessarily
          // contiguous due to empty bins between them.
          if (modible_state.count % 2 == 0 && total == modible_state.count / 2) {
            modible_state.stage_2a(counter) = 1;
            int second_index = r_counter + 2;
            while (counts(second_index, counter) == 0)
              second_index++;
            modible_state.min(counter) += r_counter * increment;
            modible_state.max(counter) = modible_state.min(counter)
                + (second_index - r_counter) * increment;
            modible_state.max_2a(counter) = modible_state.min(counter)
                + (r_counter + 1) * increment;
            modible_state.min_2a(counter) = modible_state.min(counter)
                + (second_index - 1) * increment;
            median_counts(counter) =
                accu(counts(span(r_counter + 1, second_index), counter));

          // Otherwise, proceed in an obvious fashion, i.e. keep the one bin that
          // contains the median or the two entries needed to compute it.
          } else {
            modible_state.min(counter) += r_counter * increment;
            modible_state.max(counter) = modible_state.min(counter) + increment;
            median_counts(counter) = counts(r_counter + 1, counter);
          }
          modible_state.left_counts(counter) =
              accu(counts(span(0, r_counter), counter));
        }
      }
      modible_state.range = modible_state.max - modible_state.min;
      modible_state.final_iteration = all(median_counts < <?=$sortThreshold?>);
      modible_state.iteration++;
      return true;

    } else {
      for (int counter = 0; counter < <?=$dimension?>; counter++) {
        if (modible_state.stage_2a(counter)) {
          medians(counter) = 0.5
              * (modible_state.min_2a(counter) + modible_state.max_2a(counter));
        } else {
          count = modible_state.count;
          rowvec sorted = sort(data(counter, span(0, bin_counts(counter) - 1)));
          double left_count = modible_state.left_counts(counter);
          medians(counter) = modible_state.count % 2 == 1
            ? sorted((count - 1) / 2 - left_count)
            : (sorted(count / 2 - left_count)
                + sorted(count / 2 - 1 - left_count)) / 2;
        }
      }
      modible_state.iteration++;
      return false;
    }
  }

  // The use of "count == 0" is a bit convoluted.

  //In the first stage, it is non  zero as it kept track of the count of data.

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
