// http://data.princeton.edu/wws509/notes/a2.pdf
<?
require_once "grokit_base.php";

function GLM_Predict_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className          = $t_args['className'];
    $numberCoefficients = $t_args['numberCoefficients'];
    $states             = $t_args['states'];
    $gla1               = array_keys($states)[0];
?>

using namespace arma;

class <?=$className?>ConstantState {
 private:
  // Vector of coefficients which are to be predicted. The first element is the
  // intercept if there is one.
  vec::fixed<<?=$numberCoefficients?>> beta;

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState(<?=const_typed_ref_args($states)?>) {
    beta = <?=$gla1?>.GetCoefficients();
  }

  <?=$className?>ConstantState(Json::Value state) {
    for (int counter = 0; counter < state["coefficients"].size(); counter++) {
      beta[counter] = state["coefficients"][counter].asInt();
    }
  }
};
<?
    return array(
        'kind' => 'RESOURCE',
        'name' => $className . 'ConstantState',
        'system_headers' => array('armadillo'),
        'user_headers' => array(),
        'configurable' => (count($states) == 0),
    );
}

function GLM_Predict(array $t_args, array $inputs, array $outputs, array $states)
{
    // Setting output type
    array_set_index($outputs, 0, lookupType("base::DOUBLE"));

    // Class name is randomly generated.
    $className = generate_name("GLM_Predict");

    // Initialization of local variables from template arguments

    // This is string that is either a variable name or an expression in C++;
    $response = "({$t_args['response']})";

    // Basic conversion to C++ code
    $family         = strval($t_args['family']);
    $linkFunction   = '_' . strval($t_args['link']);
    $linkDerivative = "deriv$linkFunction";
    $linkInverse    = "inverse$linkFunction";
    $variance       = "var_$family";
    $deviance       = "deviance_$family";
    $weights        = $t_args['weights'] ? $t_args['weights'] : '1.0';
    $offset         = $t_args['offset'] ? ' + ' . $t_args['offset'] : ' ';

    // The maximum allowed change allowed for convergence; It should be 0 for
    // convergence only when the clusters are stationary.
    $epsilon        = $t_args['epsilon'];

    // The mode of convergence, absolute or relative.
    $convergence = $t_args['convergence'];
    $relativeConvergence = $convergence == "relative";

    // Decides which ShouldIterate function to use.
    $isGuassian = ($family == 'gaussian');

    // The maximum iteration for the algorithm after which it will conclude
    // regardless of convergence with a warning.
    // TODO: Make sure the warning happens.
    $maxIteration   = $t_args['maxit'];

    // Variables for use in PHP only
    $debug      = $t_args['debug'];
    $order      = $t_args['order'];
    $predictors = $t_args['predictors'];

    // Construction of C++ code for first iteration
    if ($missingStart = !$t_args['start']) {
        if ($missingEta = !$t_args['eta.start']) {
            $initialEta = "$linkFunction(mu)";
            $initialMu  = !$t_args['mu.start']
                ? $t_args['family'] == "binomial"
                    ? "(($response) * $weights + 0.5) / ($weights + 1)"
                    : 'y'
                : $t_args['mu.start'];
        } else {
            $initialEta = $t_args['eta.start'];
            $initialMu  = "$linkInverse(eta)";
        }
    }

    // Interpretation of $inputs and error checks.

    // Valid family and link checking
    // An array of proper link function matched to their respective families.
    // Not all of these names are the same names as R families/links.
    // Preprocessing in rgrokit takes care of conversion.
    $properLinks = array(
         'gaussian'        => array('identity'),
         'poisson'         => array('identity', 'log',   'sqrt'),
         'gamma'           => array('identity', 'log',   'inverse'),
         'binomial'        => array('cloglog',  'logit', 'probit', 'cauchit'),
         'inverseGaussian' => array('inverseSquared')
    );

    // This checks that the family is one of the above keys.
    grokit_assert(array_key_exists($family, $properLinks),
        "Improper family: $family");

    // This checks that the array associated with the family contains the link.
    grokit_assert(in_array(substr($linkFunction, 1), $properLinks[$family]),
        "Improper link function: $linkFunction");

    // This constructs the PHP arrays necessary to process the model.
    $terms = Construct_Input_Terms($inputs, $predictors, $order);
    Remove_Redundant_Keys($terms);

    // Whether or not the model will contain an intercept.
    $interceptNeeded = $t_args['intercept'];

    // If any term contains only factors, then the intercept is not needed.
    for($counter = 0; $counter < count($terms) && $interceptNeeded; $counter ++)
        if(sizeof($terms[$counter][0]) == 0)
            $interceptNeeded = false;

    // This finishes processing the models.
    $inputVector = writeInputVector($terms, $interceptNeeded);

    // The number of coefficients to predict
    $numberCoefficients = $inputVector[0];

    // To individually set the Armadillo vector elements
    $setterCode = $inputVector[1];

    // Construction of C++ code to assign the initial coefficients if they were
    // specified. Otherwise an empty string will be produced.
    $initialCoefficientsCode = $missingStart
        ? ''
        : '{' . implode(', ', $t_args['start']);
?>

using namespace arma;

class <?=$className?>;

<?  $constantState = lookupResource(
        "statistics::GLM_Predict_Constant_State",
        array(
            'states'             => $states,
            'className'          => $className,
            'numberCoefficients' => $numberCoefficients,
        )
    );
?>


class <?=$className?> {
 private:
  // The typical constant state for an iterable GLA.
  const <?=$constantState?> & constant_state;

  // The following variables are used within the rank-1 updates. Their names
  // correspond to those found in Princeton paper for ease of comparison.
  // They are member variables rather than local variables in order to avoid
  // repeatedly allocating the same memory for them.

  // A vector that contains the current data and whose elements are individually
  // set at the beginning of the AddItem call.
  vec::fixed<<?=$numberCoefficients?>> x;

  // The response of the model. It is either a simple variable matching
  // a parameter of AddItem or an expression containing multiple parameters.
  // The details are taken care in the PHP cogeneration, which allows this to be
  // specified as any proper c++ mathematical expression.
  // TODO: Feed this through something to add namespaces, possibly RToCPP
  double y;

  // The linear predictor. It is the dot product of the coefficients and the
  // parameters of AddItem. See IBM's algorithm for more details about how it
  // is computed during the first iteration before the coefficients are first
  // predicted.
  double eta;

  // The predicted mean. This is usually found as the inverse link of eta.
  double mu;

 public:
  <?=$className?>(const <?=$constantState?>& state)
      : constant_state(state),
        x() {
  }

  bool ProcessTuple(<?=process_tuple_args($inputs, $outputs)?>) {
    // Processing of data
    x.fill(0.0);
    <?=$setterCode, PHP_EOL?>

    eta = dot(x, constant_state.beta)<?=$offset?>;
    mu = <?=$linkInverse?>(eta);
    <?=array_keys($outputs)[0]?> = mu;
    return true;
  }
};

<?php
    return array(
        'kind'            => 'GT',
        'name'            => $className,
        'system_headers'  => ['math.h',
                              'armadillo',
                              'iostream',
                              'boost/math/constants/constants.hpp',
                              'boost/math/special_functions/erf.hpp',],
	'lib_headers'     => ['glmFamilies', 'ArmaJson'],
        'user_headers'    => ['json.h'],
        'iterable'        => TRUE,
        'input'           => $inputs,
        'output'          => $outputs,
        'result_type'     => 'single',
        'generated_state' => $constantState,
        'libraries'       => ['armadillo'],
    );
}
?>
