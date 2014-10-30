<?php
// TODO: Replace more inline php with constant C++ values like RSV sampling.
// TODO: Implement triangle ineq speed up. This won't work for non-metric p < 1. Nvm, This doesn't seem to help.
// http://cseweb.ucsd.edu/~elkan/kmeansicml03.pdf
require_once "grokit_base.php";

function K_Means_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className          = $t_args['className'];
    $numberClusters     = $t_args['numberClusters'];
    $dimension          = $t_args['dimension'];
    $initialCentersCode = $t_args['initialCentersCode'];
?>
using namespace arma;
using namespace std;

class <?=$className?>ConstantState{
 private:
  // The current iteration which is to be compared to $maxIteration.
  long iteration;
  // A vector of armadillo vecs that represent the center of each cluster.
  vector<vec> centers;

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState()
      : iteration(0),
        centers(<?=$initialCentersCode?>) {
  }
};
<?php
    return array(
        'kind' => 'RESOURCE',
        'name' => $className . 'ConstantState',
        'system_headers' => array('armadillo', 'vector'),
        'user_headers' => array(),
    );
}

function K_Means(array $t_args, array $inputs, array $outputs)
{
    // Class name randomly generated
    $className = generate_name("KMC");

    // Initialization of local variables from template arguments

    $numberClusters = $t_args['number.clusters'];

    // The parameter p for the Minkowski distance metric used.
    $pMinkowski = $t_args['p'];

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

    // The following is used to choose how initial clusters are formed.

    // The sample size is used to choose initial centers. This is only necessary
    // as a template argument if kmeans++ was chosen. A warning is thrown
    // beforehand if it needlessly specified.
    if ($kMeansPP = ($t_args['centers'] == 'kMeansPP')) {
        $sampleSize = $t_args['sample.size'];
        $initialCentersCode = "$numberClusters, vec($dimension)";
    // Or a random sample of size $numberClusters is used to pick centers.
    } else if ($standardStart = ($t_args['centers'] == 'standard')) {
        $sampleSize = $numberClusters;
        $initialCentersCode = "$numberClusters, vec($dimension)";
    // Or the centers are explicitly selected.
    } else {
        grokit_assert(is_array($t_args['centers']),
            "Incorrect specification for centers.");
        $initialCentersCode = is_array($t_args['centers'][0])
            ? '{' . eval_implode('"{" . implode(", ", $val1) . "}"',
                ', ' . PHP_EOL, $t_args['centers']) . '}'
            : '{' . eval_implode('"{" . $val1 . "}"',
                ', ' . PHP_EOL, $t_args['centers']) . '}';
    }

    $codeArray = array_combine(array_keys($inputs), range(0, $dimension - 1));

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

class <?=$className?>Cluster {
 private:
  // The total number of elements in the cluster.
  // Incremented by the AddItem class function.
  long count;

  // The current center of this cluster.
  vec center;

  // The sum of the element in this cluster; used to calculate the centroid.
  vec sum;

 public:
  friend class <?=$className?>;

  <?=$className?>Cluster(const vec & _center)
      : count(0),
        center(_center),
        sum(<?=$dimension?>) {
    sum.zeros();
  }

  // This function recalculates the statistics of the cluster upon adding a
  // point, i.e. the data item. It is call AddItem only for sake of consistency
  // and this name is not relevant to Datapath - this function is only ever
  // called by AddItem of the GLA below.
  void AddItem(vec point) {
    count ++;
    sum += point;
  }

  // Combines the statistics for two clusters. These clusters should not have
  // two different centers; rather, they should be two equivalent clusters from
  // two different states of the GLA below. Like AddItem, the name is for the
  // sake of consistency and is not relevant to Datapath - this function is only
  // ever called by AddState of the GLA below.
  void AddState(<?=$className?>Cluster & other) {
    count += other.count;
    sum   += other.sum;
  }

  // Calculation of centroid in order to reposition the center of this cluster
  // for the next iteration.
  vec GetResult() {
    return sum / count;
  }

  // Distance from a given point to the center of the cluster. Used to decide
  // which cluster to move a point to in the iterative step. Based on the
  // armadillo implementation of norm which is pre-built to calculate Minkowski
  // distance.
  double Distance(vec point) {
    return norm(center - point, <?=$pMinkowski?>);
  }
};


class <?=$className?>;

<?  $constantState = lookupResource(
        "statistics::K_Means_Constant_State",
        array(
            'className'          => $className,
            'dimension'          => $dimension,
            'initialCentersCode' => $initialCentersCode,
            'numberClusters'     => $numberClusters,
        )
    ); ?>

class <?=$className?> {
 private:
  // The typical constant state for an iterable GLA.
  const <?=$constantState?> & constant_state;

  // A vector of clusters (defined above). This is contained in a vector
  // as it is unlikely that the number of clusters will be overly large.
  vector<<?=$className?>Cluster> clusters;

  // The sum of the distances between each point and its nearest cluster center.
  // This is a measure of how good the fit of the clusters is.
  // TODO: Implement a better measure, such as explained variance.
  double total_score;

  // The number of items processed so far. It is not used except in the output
  // to count the number of items processed.
  int count;

  // The counter for the GetNextResult function; denotes how many items have
  // been returned so far and is compared to $numberClusters.
  double return_counter;

  // A vector whose elements will be individually set to contain the current
  // item's data. The data is packaged as a vector due to the reliance on
  // armadillo. This exists as a member variable as opposed to a local variable
  // in order to avoid repeated memory allocation.
  vec item;

  uword index;
<?  if ($startsNeeded = ($kMeansPP || $standardStart)) {
        if ($kMeansPP) { ?>

  // The reservoir sampling GLA that is used to construct the sample which
  // k-means++ will be performed on. k-means++ is performed on a sample because
  // it requires a lot of passes over the data, 2 * ($numberClusters - 1).
<?      } else { ?>

  // The reservoir sampling GLA that is used to construct the sample which will
  // be used as the initial cluster centers.
<?      } ?>
  <?=$samplingClass?> sampling_GLA;
<?  } ?>

 public:
  <?=$className?>(const <?=$constantState?> & state)
      : constant_state(state),
        total_score(0),
        count(0),
        index(0),
        item((uword)<?=$dimension?>),
        return_counter(0) {
<?  if($startsNeeded) { ?>
    if (constant_state.iteration)
      for (int counter; counter < <?=$numberClusters?>; counter ++)
        clusters.push_back(
            <?=$className?>Cluster(constant_state.centers[counter])
        );
<?  } else { ?>
    for (int counter; counter < <?=$numberClusters?>; counter ++)
      clusters.push_back(
          <?=$className?>Cluster(constant_state.centers[counter])
      );
<?  } ?>
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    // TODO: Add variance to measure fit
<?  if($startsNeeded) { ?>
    if (!constant_state.iteration) {
      sampling_GLA.AddItem(<?=args($inputs)?>);
      return;
    }
<?  } ?>
<? foreach ($codeArray as $name => $counter) { ?>
    item[<?=$counter?>] = <?=$name?>;
<?  } ?>
    int min_cluster = 0;
    double min_distance = clusters[0].Distance(item);
    for (int counter = 1; counter < <?=$numberClusters?>; counter ++) {
        double distance = clusters[counter].Distance(item);
        if (distance < min_distance) {
          min_cluster = counter;
          min_distance = distance;
        }
    }
    clusters[min_cluster].AddItem(item);
    total_score += min_distance;
    count ++;
  }

  void AddState(<?=$className?> & other) {
<?  if ($startsNeeded) { ?>
    // The other GLA must have the same state as this one, so there is no need
    // to check if it has starts.
    if (constant_state.iteration == 0) {
      sampling_GLA.AddState(other.sampling_GLA);
      return;
    }
<?  } ?>
    for (int counter = 0; counter < <?=$numberClusters?>; counter ++)
      clusters[counter].AddState(other.clusters[counter]);
  }

  bool ShouldIterate(<?=$constantState?> & modible_state) {
<?  if ($startsNeeded) { ?>
    if (constant_state.iteration == 0) {
      // Processing of the sample into vectors. This is done here, rather than
      // in the sampling GLA, because it is possible for not only doubles to be
      // used for k-means. Therefore, it would not make sense to make a version
      // of the sample GLA use vectors instead of tuples.
      // TODO: Implement factors for k-means.

      vector<vec> sample(<?=$sampleSize?>, vec(<?=$dimension?>));
      for (int counter = 0; counter < <?=$sampleSize?>; counter ++) {
<?      foreach (range(0, $dimension - 1) as $counter) { ?>
        sample[counter][<?=$counter?>]
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
      // The total amount of initial clusters generated so far.
      int num_clusters = 1;
      // The first cluster for k-means++ is chosen at random.
      long first_index = random_index(generator);
      clusters.push_back(<?=$className?>Cluster(sample[first_index]));
      cout << sample[first_index] << endl;
      modible_state.centers[0] = sample[first_index];
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

      for (num_clusters; num_clusters < <?=$numberClusters?>; num_clusters ++) {
        total_distance = 0;
        // TODO: Improve this part, maybe with a hash table
        for (s_counter = 0; s_counter < <?=$sampleSize?>; s_counter ++) {
          min_distance = clusters[0].Distance(sample[s_counter]);
          for (c_counter = 1; c_counter < num_clusters; c_counter ++) {
            double distance = clusters[c_counter].Distance(sample[s_counter]);
            if (distance < min_distance) {
              min_cluster = c_counter;
              min_distance = distance;
            }
          }
          total_distance += min_distance;
          min_distances[s_counter] = min_distance;
        }
        random_distance = generate_canonical<double, 40>(generator) * total_distance;
        for (--s_counter; s_counter >= 0; s_counter--) {
          total_distance -= min_distances[s_counter];
          if (total_distance <= random_distance)
            break;
        }
        modible_state.centers[num_clusters] = sample[s_counter];
        clusters.push_back(<?=$className?>Cluster(sample[s_counter]));
      }
      modible_state.iteration++;
      return true;
<?      } else if ($standardStart) { ?>
      for (int counter = 0; counter < <?=$sampleSize?>; counter ++)
        modible_state.centers[counter] = sample[counter];
      modible_state.iteration++;
      return true;
<?      } ?>
    }
<?  } ?>
    bool has_converged = true;
    vec new_center(<?=$dimension?>);
    for (int counter = 0; counter < <?=$numberClusters?>; counter++) {
      new_center = clusters[counter].sum / clusters[counter].count;
      has_converged  = has_converged
                    && HasConverged(modible_state.centers[counter], new_center);
      modible_state.centers[counter] = new_center;
    }
    modible_state.iteration ++;
    return    !has_converged
           && modible_state.iteration < <?=$maxIteration?>;
  }

  bool HasConverged(vec old_center, vec new_center) {
<?  if ($relativeConvergence) { ?>
    return max(abs((old_center - new_center) / old_center)) < <?=$epsilon?>;
<?  } else { ?>
    return max(abs(old_center - new_center)) <= <?=$epsilon?>;
<?  } ?>
  }

  // TODO: Implement JSON return or something, not just printing
  void GetResult(<?=typed_ref_args($outputs)?>) {
    Json::Value result(Json::objectValue);
    result["iteration"] = (Json::Value::Int64) constant_state.iteration;
    for (int counter = 0; counter < <?=$numberClusters?>; counter ++) {
<?  foreach ($codeArray as $name => $counter) { ?>
      result["coefficients"]["<?=$name?>"].append(
        constant_state.centers[counter][<?=$counter?>]);
<?  }  ?>
    }
    cout << result << endl;
  }
};

<?
    return array(
        'kind'             => 'GLA',
        'name'             => $className,
        'system_headers'   => array('armadillo', 'vector', 'string', 'iostream'),
        'user_headers'     => array(),
        'iterable'         => TRUE,
        'input'            => $inputs,
        'output'           => $outputs,
        'result_type'      => 'single',
        'generated_states' => array($constantState),
    );
}
?>
