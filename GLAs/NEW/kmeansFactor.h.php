<?php
require_once "grokit_base.php";

function K_Means_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className          = $t_args['className'];
    $numClusters        = $t_args['numClusters'];
    $dimension          = $t_args['dimension'];
    $initialCentersCode = $t_args['initialCentersCode'];
    $numNumeric         = $t_args['numNumeric'];
    $numFactors         = $t_args['numFactors'];
    $maxPower           = $t_args['maxPower'];
?>
using namespace arma;
using namespace std;

class <?=$className?>ConstantState{
 private:
  // The current iteration which is to be compared to $maxIteration.
  long iteration;

  // Each column is the center of a cluster.
  mat::fixed<<?=$numNumeric?>, <?=$numClusters?>> centers;

  // This contains information regarding the distribution of factors in the
  // data. Rather than just containing the mode, it is organized as:
  // 1. Each slice corresponds to a factor.
  // 2. Each column corresponds to a cluster.
  // 3. Each element is 1 - f, where f is the frequency of that level of the
  //    corresponding factor for that cluster.
  cube::fixed<<?=$maxPower?>, <?=$numClusters?>, <?=$numFactors?>> modes;
<?  if ($normalized) ?>

  // This contains 3 column - max, min, range - with an entry for each input.
  // Only range and minimum are needed for the linear transformation. Max is
  // stored in case of future additions.
  mat::fixed<<?=$numNumerics?>, 3> transformations;

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState()
      : iteration(0),
        centers(<?=$initialCentersCode?>) {
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

function K_Means(array $t_args, array $inputs, array $outputs)
{
    // Class name randomly generated
    $className = generate_name("KMC");

    // Initialization of local variables from template arguments

    $numClusters   = $t_args['number.clusters'];
    $debug         = $t_args['debug'];
    $normalization = $t_args['normalization'];

    // The parameter p for the Minkowski distance metric used.
    $pMinkowski = $t_args['p'];
    grokit_assert($pMinkowski != 0, "The parameter p cannot be exactly zero.");

    // The parameter m for the fuzzified membership used.
    $mFuzzifier = $t_args['m'];
    grokit_assert($pMinkowski >= 1, "The parameter m must be at least one.");

    // The dimension of the data, i.e. how many elements are in each item.
    $dimension = count($inputs);

    // The max iteration for the algorithm after which it will conclude
    // regardless of convergence with a warning.
    // TODO: Make sure the warning happens.
    $maxIteration = $t_args['max.iteration'];

    // The maximum allowed change allowed for convergence; It should be 0 for
    // convergence only when the clusters are stationary.
    $epsilon = $t_args['epsilon'];

    // The mode of convergence, absolute or relative.
    $convergence = $t_args['convergence'];
    $relativeConvergence = $convergence == "relative";
    grokit_error_if($epsilon < 0, "Epsilon cannot be negative.");
    grokit_error_if($epsilon > 1 && $relativeConvergence,
        "Epsilon must be between 0 and 1 for relative convergence.");

    // The following is used to choose how initial clusters are formed.

    // The sample size is used to choose initial centers. This is only necessary
    // as a template argument if kmeans++ was chosen. A warning is thrown
    // beforehand if it needlessly specified.
    if ($kMeansPP = ($t_args['centers'] == 'kMeansPP')) {
        $sampleSize = $t_args['sample.size'];
        $initialCentersCode = '';
    // Or a random sample of size $numClusters is used to pick centers.
    } else if ($standardStart = ($t_args['centers'] == 'standard')) {
        $sampleSize = $numClusters;
        $initialCentersCode = '';
    // Or the centers are explicitly selected.
    } else {
        grokit_assert(is_array($t_args['centers']),
            "Incorrect specification for centers.");
        $initialCentersCode = is_array($t_args['centers'][0])
            ? '{' . eval_implode('implode(", ", $val1)',
                ', ' . PHP_EOL, $t_args['centers']) . '}'
            : '{' . implode(',' . PHP_EOL, $t_args['centers']) . '}';
    }

    $numeric = array();
    $factors = array();
    $powers  = array();
    foreach ($inputs as $name => $type)
      if ($type->is("numeric")) {
        $numeric[] = $name;
      } else {
        $factors[] = $name;
        $powers[] = $type->get("cardinality");
      }
    $numNumeric = count($numeric);
    $numFactors = count($factors);
    $maxPower   = max($powers);
    grokit_error_assert($relativeConvergence || $numFactors == 0,
                    "Relative convergence is required for categorical data.");

    $codeArray1 = array_combine(array_keys($inputs), range(0, $dimension - 1));
    $codeArray2 = array_combine($numeric, range(0, $numNumeric - 1));
    $codeArray3 = array_combine($factors, range(0, $numFactors - 1));

    if ($kMeansPP || $standardStart) {
        // The GLA for reservoir sampling is created.
        $samplingClass = lookupGLA(
            "statistics::Reservoir_Sampling",
            array('threshold.coefficient' => 22,
                  'sample.size'           => $sampleSize),
            $inputs,
            $inputs
        );
    }
?>

using namespace arma;
using namespace std;

class <?=$className?>;

<?  $constantState = lookupResource(
        "statistics::K_Means_Constant_State",
        array(
            'className'          => $className,
            'dimension'          => $dimension,
            'initialCentersCode' => $initialCentersCode,
            'numClusters'     => $numClusters,
            'numNumeric'         => $numNumeric,
            'numFactors'         => $numFactors,
            'maxPower'           => $maxPower,
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
  vec::fixed<<?=$numNumeric?>> item;

  // A matrix where the k-th column is the sum of the set of points in the k-th
  // cluster. These sums are used to calculate the centroids during iteration.
  mat::fixed<<?=$numNumeric?>, <?=$numClusters?>> sums;

  // A cube where the k-th slice is the frequency of that factor for each level
  // and cluster center, where columns correspond to clusters.
  cube::fixed<<?=$maxPower?>, <?=$numClusters?>, <?=$numFactors?>> frequencies;

  // A rowvec containing the counts for each cluster. For fuzzy clustering, this
  // is not the number of items in the cluster, but the sum of the weights for
  // that cluster.
  rowvec::fixed<<?=$numClusters?>> counts;

<?  if ($normalized) { ?>
  // This contains 3 column - max, min, range - with an entry for each input.
  // Only range and minimum are needed for the linear transformation. Max is
  // stored in case of future additions.
  mat::fixed<<?=$numNumerics?>, 3> transformations;

<?  } ?>
  // The sum of the distances between each point and its nearest cluster center.
  // This is a measure of how good the fit of the clusters is.
  // TODO: Implement a better measure, such as explained variance.
  double total_score;

<?  if ($startsNeeded = ($kMeansPP || $standardStart)) {
        if ($kMeansPP) { ?>
  // The reservoir sampling GLA that is used to construct the sample which
  // k-means++ will be performed on. k-means++ is performed on a sample because
  // it requires a lot of passes over the data, 2 * ($numClusters - 1).
<?      } else { ?>
  // The reservoir sampling GLA that is used to construct the sample which will
  // be used as the initial cluster centers.
<?      } ?>
  <?=$samplingClass?> sampling_GLA;

  // This is used to initialize the frequencies for the factors.
  // The columns correspond to the levels of a factor, with a column per factor.
  mat::fixed<<?=$maxPower?>, <?=$numFactors?>> overall_frequencies;

<?  } ?>
  // These four variables are used in AddItem to fully utilize Armadillo:

<?  if ($mFuzzifier > 1) { ?>
  // This is used to store the denominator of the sum used for fuzzy clustering.
  // It is the sum of the distances between each center and the current point,
  // each raised to the p / (m - 1) power.
  double distance_sum;
<?  } else { ?>
  // The index of the cluster center which is closest to the point.
  uword index;
<?  } ?>

  // Stores the distance from each center to the current point, the sum of
  // shifted centers. The 1/p root isn't taken but this is does not affect the
  // computation of the index of the closest points. It is further used in fuzzy
  // clustering for computing the weights, if applicable.
  rowvec::fixed<<?=$numClusters?>> distances;

  // Stores each center in a column after subtracting the current point from
  // each column and later the various powers of itself.
  mat::fixed<<?=$numNumeric?>, <?=$numClusters?>> shifted_centers;

  // Stores the distances for the factors. Acts much like shifted centers.
  mat::fixed<<?=$numFactors?>, <?=$numClusters?>> shifted_modes;

 public:
  <?=$className?>(const <?=$constantState?> & state)
      : constant_state(state),
        item(),
        sums(zeros<mat>(<?=$dimension?>, <?=$numClusters?>)),
        counts(zeros<rowvec>(<?=$numClusters?>)),
<?  if ($normalized) { ?>
        transformations()
<?  } ?>
        total_score(0),
<?  if ($startsNeeded) { ?>
        sampling_GLA(),
        overall_frequencies(zeros<mat>(<?=$maxPower?>, <?=$numFactors?>)),
<?  } ?>
<?  if ($mFuzzifier > 1) { ?>
        distance_sum(0),
<?  } else { ?>
        index(0),
<?  } ?>
        distances(),
        shifted_centers() {
<?  if ($normalized) { ?>
    transformations.col(0).fill(-numeric_limits<double>::infinity());
    transformations.col(1).fill( numeric_limits<double>::infinity());
<?  } ?>
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    // TODO: Add diagonistic statistics to see how good the fit is.
<?  if ($normalized) { ?>
    if (constant_state.iteration == 0) {
<? foreach ($codeArray2 as $name => $counter) { ?>
      if (<?=$name?> > extremities(<?=$counter?>, 0))
        extremities(<?=$counter?>, 0) = <?=$name?>;
      if (<?=$name?> < extremities(<?=$counter?>, 1))
        extremities(<?=$counter?>, 1) = <?=$name?>;
<?  } ?>
    }

<?  } ?>
<?  if($startsNeeded) { ?>
    if (constant_state.iteration == 0) {
      sampling_GLA.AddItem(<?=args($numeric)?>);
<?  foreach ($codeArray2 as $name => $counter) { ?>
      overall_frequencies[<?=$name?>][<?=$counter?>]++;
<?  } ?>
      return;
    }

<?  } ?>
<? foreach ($codeArray2 as $name => $counter) { ?>
    item[<?=$counter?>] = <?=$name?>;
<?  } ?>
<?  if ($normalized) { ?>
    item -= transformations.col(1);
    item /= transformations.col(2);
<?  } ?>
    shifted_centers = constant_state.centers;
    shifted_centers.each_col() -= item;
<?  foreach ($codeArray2 as $factor => $counter) { ?>
    shifted_modes.row(<?=$counter?>) = constant_state.frequencies
        .slice(<?=$counter?>).row(<?=$factor?>);
<?  } ?>
<?  if ($pMinkowski % 2 != 0) { ?>
    shifted_centers = abs(shifted_centers);
<?  } ?>
<?  if ($pMinkowski % 1 == 0 && -8 <= $pMinkowski && $pMinkowski <= 10) { ?>
<?      if ($pMinkowski < 0) { ?>
    shifted_centers = 1 / shifted_centers;
    shifted_modes = 1 / shifted_modes;
<?      }
        for ($counter = 2; $counter <= $pMinkowski; $counter ++) { ?>
    shifted_centers %= shifted_centers;
    shifted_modes %= shifted_modes;
<?      } ?>
<?      for ($counter = -2; $counter >= $pMinkowski; $counter --) { ?>
    shifted_centers /= shifted_centers;
    shifted_modes /= shifted_modes;
<?      } ?>
<?  } else { ?>
    shifted_centers = pow(shifted_centers, <?=$pMinkowski?>);
    shifted_modes = pow(shifted_modes, <?=$pMinkowski?>);
<? } ?>
    distances  = sum(shifted_centers, 0);
    distances += sum(shifted_modes, 0);

<?  if ($mFuzzifier == 1) {
        if ($pMinkowski < 0) { ?>
    total_score += distances.max(index);
<?      } else { ?>
    total_score += distances.min(index);
<?      } ?>
    sums.col(index) += item;
    counts[index] ++;
<?      foreach ($codeArray2 as $factor => $counter) { ?>
    frequencies(<?=$factor?>, index, <?=$counter?>) ++;
<?      } ?>
<?      if ($debug) { ?>
    if (sum(counts) <= 10) {
      cout << "Item number: " << sum(counts) << endl;
      cout << "Item: " << endl << item << endl;
      cout << "Added to cluster " << index << endl;
      cout << "Current counts:" << endl << counts << endl;
      cout << "Current sums:" << endl << sums << endl;
      cout << "Current frequencies:" << endl << frequencies << endl;
    }
<?      }
    } else { ?>
    distances = pow(distances + 0.000000001, <?= -1 / ($mFuzzifier - 1) ?>);
    distances /= sum(distances);
    counts += distances;
    sums   += item * distances;
<?      foreach ($codeArray2 as $factor => $counter) { ?>
    frequencies.slice(<?=$counter?>).row(<?=$factor?>) += distances;
<?      } ?>
<?      if ($debug) { ?>
    if (sum(counts) <= 10) {
      cout << "Item number: " << sum(counts) << endl;
      cout << "Item:" << endl << item << endl;
      cout << "Weights:" << endl << distances << endl;
      cout << "Current counts:" << endl << counts << endl;
      cout << "Current sums:" << endl << sums << endl;
      cout << "Current frequencies:" << endl << frequencies << endl;
    }
<?      }
    } ?>
  }

  void AddState(<?=$className?> & other) {
<?  if ($startsNeeded) { ?>
    // The other GLA must have the same state as this one, so there is no need
    // to check if it has starts.
    if (constant_state.iteration == 0) {
      sampling_GLA.AddState(other.sampling_GLA);
      overall_frequencies += other.overall_frequencies;
      return;
    }
<?  }
    if ($debug) { ?>
    cout << "Before Sums:" << endl << sums << endl;
    cout << "Before Counts:" << endl << counts << endl;
    cout << "Before Frequencies:" << endl << frequencies << endl;
<?  } ?>
    sums   += other.sums;
    counts += other.counts;
    frequencies += other.frequencies
<?  if ($debug) { ?>
    cout << "After Sums:" << endl << sums << endl;
    cout << "After Counts:" << endl << counts << endl;
    cout << "After Frequencies:" << endl << frequencies << endl;
<?  } ?>
  }

  bool ShouldIterate(<?=$constantState?> & modible_state) {
<?  if ($debug) { ?>
    cout << "Iteration: " << modible_state.iteration + 1 << endl;
    cout << "Sums:" << endl << sums << endl;
    cout << "Counts:" << endl << counts << endl;
    cout << "Frequencies:" << end << frequencies << endl;
<?  } ?>
<?  if ($normalized) { ?>
    if (constant_state.iteration == 0) {
      transformations.col(2) = transformations.col(0) - transformations.col(1);
      modible_state.transformations = transformations;
    }
<?  } ?>
<?  if ($startsNeeded) { ?>
    if (constant_state.iteration == 0) {
      // Processing of the sample into vectors. This is done here, rather than
      // in the sampling GLA, because it is possible for not only doubles to be
      // used for k-means. Therefore, it would not make sense to make a version
      // of the sample GLA use vectors instead of tuples.
      // TODO: Implement factors for k-means.

      mat::fixed<<?=$numNumeric?>, <?=$sampleSize?>> sample;
      for (int counter = 0; counter < <?=$sampleSize?>; counter ++) {
<?      foreach ($codeArray2 as $counter) { ?>
        sample.row(<?=$counter?>)[counter]
          = get<<?=$counter?>>(sampling_GLA.GetSample(counter));
<?      }  ?>
      }
<?      if ($kMeansPP) { ?>
      // Implementation of k-means++

      // A random device used to generate the seed for the generator
      std::random_device random_seed;

      // A random engine used to generated random variates for the above.
      std::default_random_engine generator(random_seed());

      // A uniform discrete distribution that generates a random variate between
      // 0 and sample size - 1, inclusive. This is used to randomly decide which
      // element of sample should be used as the first cluster center.
      std::uniform_int_distribution<long> random_index(0, <?=$sampleSize?> - 1);

      // The first cluster for k-means++ is chosen at random.
      long first_index = random_index(generator);
      modible_state.centers.col(0) = sample.col(first_index);

      // The total amount of initial clusters generated so far.
      int num_clusters = 1;

      // These are declared here only to avoid re-allocation of memory.
      // Otherwise, they would be declared within the for loop.

      // Sum of the distances between each point and its nearest cluster center.
      int total_distance;
      // Stores the distance between each point and its nearest cluster center.
      vec min_distances(<?=$sampleSize?>);
      // The distance between the current point and its nearest cluster center.
      double min_distance;
      // The index of the closest cluster to the current point.
      long min_cluster;
      // Iterates over the sample.
      long s_counter;
      // Iterates over the cluster centers.
      long c_counter;
      // A random double between 0 and the total distance. It is used to select
      // the next cluster center.
      double random_distance;

      for (num_clusters; num_clusters < <?=$numClusters?>; num_clusters ++) {
        total_distance = 0;
        // TODO: Improve this part, maybe with a hash table
        for (s_counter = 0; s_counter < <?=$sampleSize?>; s_counter ++) {
          min_distance = norm(
                modible_state.centers.col(0) - sample(s_counter),
                <?=$pMinkowski?>);
          for (c_counter = 1; c_counter < num_clusters; c_counter ++) {
            double distance = norm(
                modible_state.centers.col(c_counter) - sample.col(s_counter),
                <?=$pMinkowski?>);
            if (distance < min_distance) {
              min_cluster = c_counter;
              min_distance = distance;
            }
          }
          total_distance += min_distance;
          min_distances[s_counter] = min_distance;
        }
        random_distance = generate_canonical<double, 40>(generator)
                        * total_distance;
        for (--s_counter; s_counter >= 0; s_counter--) {
          total_distance -= min_distances[s_counter];
          if (total_distance <= random_distance)
            break;
        }
        modible_state.centers.col(num_clusters) = sample.col(s_counter);
      }
<?      } else if ($standardStart) { ?>
      for (int counter = 0; counter < <?=$sampleSize?>; counter ++)
        modible_state.centers[counter] = sample[counter];
<?      } ?>
      overall_frequencies /= sum(overall_frequencies.col(0));
      int f_counter = 0, l_counter = 0, c_counter = 0;
      for (f_counter; f_counter < <?=$numFactors?>; f_counter++)
        for (l_counter; l_counter < <?=$numFactors?>; l_counter++)
          for (c_counter; c_counter < <?=$numFactors?>; c_counter++)
            modible_state.modes(l_counter, c_counter, f_counter) =
              overall_frequencies(l_counter, f_counter);
      modible_state.iteration++;
      return true;
    }

<?  } ?>
<?  if (!normalized) { ?>
    sums.each_col() *= modible_state.transformations.col(2);
    sums.each_col() += modible_state.transfomrations.col(1);
<?  } ?>
    bool has_converged = true;
    vec::fixed<<?=$numNumeric?>> new_center();
    vec::fixed<<?=$numNumeric?>> new_mode();

    for (int c_counter = 0; c_counter < <?=$numClusters?>; c_counter++) {
      if (counts[c_counter] > 0)
        new_center = sums.col(c_counter) / counts[c_counter];
      else
        new_center = modible_state.centers.col(c_counter);
      has_converged  = has_converged
                    && HasConverged(modible_state.centers.col(c_counter),
                                    new_center);
      modible_state.centers.col(c_counter) = new_center;
    }

    for (int f_counter = 0; f_counter < <?=$numFactors?>)
      for (c_counter = 0; c_counter < <?=$numClusters?>; c_counter ++) {
        if (counts[c_counter] > 0)
          new_mode = 1 - frequencies.slice(f_counter).col(c_counter)
                   / sum(frequencies.slice(0).col(0));
        else
          new_mode = modible_state.modes.slice(f_counter).col(c_counter);
        has_converged  = has_converged
                      && HasConverged(
                             modible_state.modes.slice(f_counter).col(c_counter),
                             new_mode);
        modible_state.modes.slice(f_counter).col(c_counter) = new_mode;
    }

    modible_state.iteration ++;
    return    !has_converged
           && modible_state.iteration < <?=$maxIteration?>;
  }

  bool HasConverged(vec old_center, vec new_center) {
<?  if ($relativeConvergence) { ?>
    return max(abs((old_center - new_center) / old_center)) <= <?=$epsilon?>;
<?  } else { ?>
    return max(abs(old_center - new_center)) <= <?=$epsilon?>;
<?  } ?>
  }

  // TODO: Implement JSON return or something, not just printing
  void GetResult(<?=typed_ref_args($outputs)?>) {
    Json::Value result(Json::objectValue);
    result["iteration"] = (Json::Value::Int64) constant_state.iteration;
    for (int counter = 0; counter < <?=$numClusters?>; counter ++) {
<?  foreach ($codeArray1 as $name => $counter) { ?>
      result["centers"][<?=$counter?>].append(
        constant_state.centers.col(counter)[<?=$counter?>]);
<?  }  ?>
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
