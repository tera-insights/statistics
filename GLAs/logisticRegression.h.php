<?
// An implementation of logistic regression using gradient descent, in constrast
// to the generalized linear model which uses Newton's Method. The trade-off is
// that gradient descent needs more iterations to converge but it requires O(d)
// space and time for d dimensions, whereas Newton's Method will require O(d^2).

// This implementation is specifically intended for use with an RGPS method that
// imitates Oracle Health Sciences'. As such, many features seen in the GLM are
// missing and the accepted input is much narrower.

// The white paper on RGPS can be found at:
// http://www.oracle.com/us/industries/health-sciences/hs-regression-adjusted-gps-wp-1949689.pdf

// The typical algorithm for the gradient descent method can be found at:
// http://qwone.com/~jason/writing/lr.pdf

require_once "grokit_base.php";

function Logistic_Regression_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className    = $t_args['className'];
    $numberDrugs  = $t_args['numberDrugs'];
    $numberModels = $t_args['numberModels'];
    $numberStrata = $t_args['numberStrata'];
?>

using namespace arma;
using namespace std;

class <?=$className?>ConstantState {
 private:
  // The current iteration which is to be compared to $maxIteration.
  long iteration;

  // A matrix of coefficients. The format of the entries are as follows:
  // a) The nth column corresponds to a model for the nth AE, which is computed
  //    somewhat indepedently of the other models.
  // b) The first x slots are reserved for group intercept terms, where x is the
  //    maximum number of groups across all models.
  // c) Not all x slots are necessarily used for each model. Hence there will be
  //    padding 0s between terms corresponding to strata groups and drugs.
  // For the constant state, only the maximum number of terms is needed, which is
  // x + number of drugs. Not necessarily all drugs are used for each model, as
  // per the ELR analysis.
  // The mathematical object is separated into two different data structures to
  // ease access and computation.
  mat::fixed<<?=$numberDrugs?>, <?=$numberModels?>> beta_drugs;
  mat beta_strata;

  // The mappings from strata (represented as an integer) to its group's index
  // for each adverse effect model, which are independent of each other.
  //int groupings [<?=$numberModels?>][<?=$numberStrata?>];
  umat::fixed<<?=$numberStrata?>, <?=$numberModels?>> groupings;

  // The number of slots reserved for intercepts in any beta_drugs column. The
  // matrix cannot be jagged, so some unneeded memory is allocated.
  int num_groups;

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState()
      : iteration(0),
        beta_strata(),
        groupings(),
        num_groups(0) {
    beta_drugs.ones();
  }
};
<?php
    return array(
        'kind' => 'RESOURCE',
        'name' => $className . 'ConstantState',
        'system_headers' => array('armadillo'),
        'user_headers' => array(),
    );
}
function Logistic_Regression(array $t_args, array $inputs, array $outputs) {
    // Class name is randomly generated
    $className = generate_name("LogReg");

    // These are hard coded checks to ensure that the input consists of:
    // 1) Sex, Age, Year categorical variables that encode the strata.
    // 2) A sparse vector indicating which drugs the patient used.
    // 3) A sparse vector indicating which adverse effects the patient reported.
    // Each row of the data corresponds to a single patient.
    // TODO: Add check

    // Names of inputs/arguments are recorded to improve readability
    $sex     = array_keys($inputs)[0];
    $age     = array_keys($inputs)[1];
    $year    = array_keys($inputs)[2];
    $drugs   = array_keys($inputs)[3];
    $effects = array_keys($inputs)[4];
    $numSex  = array_get_index($inputs, 0)->get("cardinality");
    $numAge  = array_get_index($inputs, 1)->get("cardinality");
    $numYear = array_get_index($inputs, 2)->get("cardinality");

    // Maximum number of non-intercept coefficients for a given model. Most of
    // of them will be 0 for any given model due to the prior ELR analysis.
    $numberDrugs = 4194;
    $numberModels = 15092;

    // TODO: Change it from hard-coding to something else.

    // Number of possible strata. Not all necessarily appear.
    $numberStrata = 1188;
    $epsilon = 0.1;
    $alpha = .25;
    $splineLink = TRUE;
    $relativeChange = 0.05;
?>

using namespace std;
using namespace arma;

class <?=$className?>;

<?  $constantState = lookupResource(
        "statistics::Logistic_Regression_Constant_State",
        array(
            'className'    => $className,
            'numberDrugs'  => $numberDrugs,
            'numberModels' => $numberModels,
            'numberStrata' => $numberStrata,
        )
    ); ?>

class <?=$className?> {
 private:
  // The typical constant state for an iterable GLA.
  const <?=$constantState?> & constant_state;

  // Statistics used for RGPS. All names come from the cited papers.

  // The gamma parameter for the shrinkage of the stratum counts.
  // double gamma_0;

  // The counts of various strata. Some strata are never seen in the FDA data
  // but all K strata are used in the gamma-poisson shrinkage. The strata are
  // hashed to their count based on a numeral system with multiple bases.
  uvec::fixed<<?=$numberStrata?>> N_s;

  // The number of occurrences for each adverse affect.
  uvec::fixed<<?=$numberModels?>> N_plus;

  // Statistics used for typical gradient descent.

  // The sums of y_i * x_i * g(-y_i * z_i) for each adverse effect.
  mat::fixed<<?=$numberDrugs?>, <?=$numberModels?>> XYG_drugs;
  mat XYG_strata;

  // The learning rate. This is needed for typical gradient descent and maybe
  // for RGPS.
  // TODO: Figure out update for RGPS exactly.
  double epsilon;

  int tCounter = 0; // tuple counter

 public:
  <?=$className?>(const <?=$constantState?> &state)
      : constant_state(state),
        epsilon(<?=$epsilon?>),
        XYG_strata(state.num_groups, <?=$numberModels?>),
        N_s(),
        N_plus(),
        XYG_drugs() {
    cout << "Started log reg " << endl;
    N_s.zeros();
    N_plus.zeros();
    XYG_drugs.zeros();
    XYG_strata.zeros();
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    int s = HashStrata(<?=$sex?>, <?=$age?>, <?=$year?>);
    // During the first iteration we proceed with the RGPS analysis.
    // TODO: Expand implementation of RGPS.
    if (constant_state.iteration == 0) {
      N_s[s]++;
      for (int index : <?=$effects?>.GetResult())
        N_plus[index]++;
    } else {
      tCounter ++;

      arma::colvec z = 

      // first asume there are no adverse efects
      for (int counter = 0; counter < <?=$numberModels?>; counter++) {
	int g_i = constant_state.groupings(s,counter);
        double z = constant_state.beta_strata(g_i, counter) + <?=$drugs?> * constant_state.beta_drugs.unsafe_col(counter);
        double g_z = g(z);
        XYG_strata(g_i,counter) -= g_z;
        for (int index : <?=$drugs?>.GetResult())
           XYG_drugs(index,counter) -= g_z;
       }

      // correct the counters for which there are adverse effects
      for (int counter : <?=$effects?>.GetResult()) {
        int g_i = constant_state.groupings(s,counter);
        XYG_strata(g_i,counter) += 1;
        for (int index : <?=$drugs?>.GetResult())
           XYG_drugs(index, counter) += 1;
      }      

      if (tCounter % 1000 == 0){
	cout << "."; cout.flush();
      }
    }
  }

  void AddState(const <?=$className?>& other) {
    if (constant_state.iteration == 0) {
    //      cout << "Added states for iteration 1" << endl;
      N_s += other.N_s;
      N_plus += other.N_plus;
    } else {
      XYG_drugs += other.XYG_drugs;
      XYG_strata += other.XYG_strata;
    }
  }

  bool ShouldIterate(<?=$constantState?>& modible_state) {
    if (constant_state.iteration == 0) {
      cout << "ShouldIterate iteration 1" << endl;
      // TODO: Actually sort the strata before grouping by RGPS.
      int number_strata = sum(N_s);
      int max_group_counter = 0;
      // Creating the mappings from strata to group for each model.
      // TODO: Use correct ordering.
      for (int model = 0; model < <?=$numberModels?>; model++) {
        int group_count = 0;
        int group_counter = 0;
        int N_min = (int) fmax(10, sqrt(N_plus[model]) / 2);
        for (int strata = 0; strata < <?=$numberStrata?>; strata++) {
          group_count += N_s[strata];
          modible_state.groupings(strata,model) = group_counter;
          if (group_counter >= N_min && strata != <?=$numberStrata - 1?>)
            group_counter++;
        }
        if (group_counter > max_group_counter)
          max_group_counter = group_counter;
      }
      int num_groups = max_group_counter + 1;
      cout << "num groups: " << num_groups << endl;
      modible_state.beta_strata = ones(num_groups, <?=$numberModels?>);
      modible_state.beta_drugs = zeros(<?=$numberDrugs?>, <?=$numberModels?>);
      modible_state.num_groups = num_groups;
      modible_state.iteration++;
      cout << "Iteration " << modible_state.iteration << " finished." << endl;
      return true;
    } else {
      double relative_change = epsilon * fmax(
          max(max(abs(XYG_drugs) / modible_state.beta_drugs)),
          max(max(abs(XYG_strata) / modible_state.beta_strata))
      );
      modible_state.beta_drugs += epsilon * XYG_drugs;
      modible_state.beta_strata += epsilon * XYG_strata;
      modible_state.iteration++;
      cout << "Change: " << relative_change << endl;
      cout << "Iteration " << modible_state.iteration << " finished." << endl;
      return relative_change > <?=$relativeChange?>;
    }
  }

  void GetResult(<?=$outputs[array_keys($outputs)[0]]?> & <?=array_keys($outputs)[0]?>) {
<?  ProduceResult(
        array_keys($outputs)[0],
        array(
            'iteration' => 'constant_state.iteration',
            'drug coefficients' => 'constant_state.beta_drugs',
            'strata coefficients' => 'constant_state.beta_strata',
        )
    ) ?>
  }

  inline int HashStrata(<?=const_typed_ref_args(array_slice($inputs, 0, 3, TRUE))?>) {
    return <?=$sex?> * <?=$numAge * $numYear?> + <?=$age?> * <?=$numYear?> + <?=$year?>;
  }

  /* Experiments with g(z)=z*z show that this is not the bottleneck, only 10% speed improvement */
  inline double g(double z) {
<? if ($splineLink) { ?>
    return (z < 0) ? <?=2 * $alpha?> / (1 + exp(<?=2 * ($alpha - 1)?> * z))
                   : <?=2 * $alpha - 1 ?> + <?=2 * (1 - $alpha)?> / (1 + exp(<?=-2 * $alpha?> * z));
<? } else { ?>
    return 1 / (1 + exp(-z));
<? } ?>
  }
};

<?
    return array(
        'kind' => 'GLA',
        'name' => $className,
        'system_headers' => array(
            'math.h',
            'armadillo',
        ),
	'lib_headers' => ['ArmaJson'],
        'user_headers' => array('json.h'),
        'iterable' => TRUE,
        'input' => $inputs,
        'output' => $outputs,
        'result_type' => 'single',
        'generated_states' => array($constantState),
        'libraries' => array('armadillo'),
    );
}

?>
