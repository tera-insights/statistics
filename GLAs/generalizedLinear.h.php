
<?
// http://data.princeton.edu/wws509/notes/a2.pdf

function GLM_Constant_State(array $t_args)
{
    // Grabbing variables from $t_args
    $className               = $t_args['className'];
    $numberCoefficients      = $t_args['numberCoefficients'];
    $initialCoefficientsCode = $t_args['initialCoefficientsCode'];
?>

using namespace arma;

class <?=$className?>ConstantState {
 private:
  // Current iteration of the fitting.
  long iteration;

  // Vector of coefficients which are to be predicted. The first element is the
  // intercept if there is one.
  vec::fixed<<?=$numberCoefficients?>> beta;

 public:
  friend class <?=$className?>;

  <?=$className?>ConstantState()
      : iteration(0),
        beta(<?=$initialCoefficientsCode?>) {
  }
};
<?
    return [
        'kind' => 'RESOURCE',
        'name' => $className . 'ConstantState',
        'system_headers' => ['armadillo'],
        'user_headers' => [],
    ];
}

function GLM(array $t_args, array $inputs, array $outputs)
{
    // Setting output type
    $outputs = ['_output' => lookupType('JSON')];

    // Class name is randomly generated.
    $className = generate_name("GLM");

    // Processing of template arguments

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
    $offset         = $t_args['offset'] ? ' + ' . $t_args['offset'] : '';

    // The maximum allowed change allowed for convergence; It should be 0 for
    // convergence only when the clusters are stationary.
    $epsilon        = $t_args['epsilon'];

    // The mode of convergence, absolute or relative.
    $convergence = $t_args['convergence'];
    $relativeConvergence = $convergence == "relative";

    // Decides which ShouldIterate function to use.
    $isGuassian = $family == 'gaussian';

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
    $properLinks = [
         'gaussian'        => array('identity'),
         'poisson'         => array('identity', 'log',   'sqrt'),
         'gamma'           => array('identity', 'log',   'inverse'),
         'binomial'        => array('cloglog',  'logit', 'probit', 'cauchit'),
         'inverseGaussian' => array('inverseSquared')
    ];

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

    // The number of coefficients to predict.
    $numberCoefficients = $inputVector[0];

    // To individually set the Armadillo vector elements.
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
        "statistics::GLM_Constant_State",
        ['className'               => $className,
         'numberCoefficients'      => $numberCoefficients,
         'initialCoefficientsCode' => $initialCoefficientsCode,]
    );
?>


class <?=$className?> {
 private:
  // The typical constant state for an iterable GLA.
  const <?=$constantState?>& constant_state;

  // This is one of the two matrices used to compute the next iteration of the
  // coefficients. It is equal to XWX, where X is the matrix containing all data
  // and W is the weights matrix. It is computed by the rank-1 update system in
  // AddItem.
  mat::fixed<<?=$numberCoefficients?>, <?=$numberCoefficients?>> XWX;

  // The other matrix used to compute the next iteration of the coefficients.
  // It is equal to XWZ, where Z is the matrix of working dependent variables.
  vec::fixed<<?=$numberCoefficients?>> XWZ;

  // The total numbers of items processed so far.
  long count;

  // The deviance recorded during the last iteration.
  double deviance;

  // The Pearson's chi squared statistic recorded during the last iteration.
  double pearson;

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

  // The working dependant variable.
  double z;

  // The weight for the current item.
  double w;

 public:
  <?=$className?>(const <?=$constantState?>& state)
      : constant_state(state),
        count(0),
        x(),
        deviance(0),
        pearson(0),
        XWX(zeros<mat>(<?=$numberCoefficients?>, <?=$numberCoefficients?>)),
        XWZ(zeros<vec>(<?=$numberCoefficients?>)) {
  }

  void AddItem(<?=const_typed_ref_args($inputs)?>) {
    // Processing of data
    ++count;
    x.fill(0.0);
    <?=$setterCode, PHP_EOL?>
    y = <?=$response?>;

    // Rank-1 updates
<?  if ($missingStart) { ?>
    // Initial values will be found by IBM's algorithm.
    if (constant_state.iteration == 0) {
<?      if ($missingEta) { ?>
      // Eta is built off mu.
      mu  = <?=$initialMu?>;
      eta = <?=$initialEta?>;
<?      } else { ?>
      // Mu is built off eta.
      eta = <?=$initialEta?>;
      mu  = <?=$initialMu?>;
<?      } ?>
    } else {
      eta = dot(x, constant_state.beta)<?=$offset?>;
      mu = <?=$linkInverse?>(eta);
    }
<?  } else { ?>
    eta = dot(x, constant_state.beta)<?=$offset?>;
    mu = <?=$linkInverse?>(eta);
<?  } ?>
    z = eta + (y - mu) * <?=$linkDerivative?>(mu);
    w = <?=$weights?>
      / (<?=$variance?>(mu) * std::pow(<?=$linkDerivative?>(mu), 2));

    XWZ += x * w     * z;  // Update of XWZ
    XWX += x * x.t() * w;  // Update of XWX

    deviance += <?=$weights?> * <?=$deviance?>(y, mu);
    pearson += <?=$weights?> * (y - mu) * (y - mu) / <?=$variance?>(mu);

<?  if ($debug) { ?>
    // Debug statements
    if (count <= 10){
      cout << "x: " << endl << x << endl;
      cout << "y: " << y << " eta: " << eta << " mu: " << mu
           << " z: " << z << " w: " << w
           << " deivance: " << deviance
           << " pearson: " << pearson << endl;
      cout << "XWZ: " << endl << XWZ << endl;
      cout << "XWX: " << endl << XWX << endl;
    };
<?  } ?>
  }

  // Straightforward since all the state is numerically additive
  void AddState (const <?=$className?>& other) {
<?  if ($debug) { ?>
    // Debug statements
    cout << "Merging" << endl;
    cout << "XWX:" << endl << XWX << endl;
    cout << "XWZ:" << endl << XWZ << endl;

<?  } ?>
    count += other.count;
    XWX   += other.XWX;
    XWZ   += other.XWZ;

    deviance += other.deviance;
    pearson += other.pearson;
<?  if ($debug) { ?>
    // Debug statement
    cout << "After Merge" << endl;
    cout << "XWX:" << endl << XWX << endl;
    cout << "XWZ:" << endl << XWZ << endl;
<?  } ?>
  }

  // The next iteration of the coeffficient vector is solved by computing the
  // generalized inverse (Moore Penrose pseudoinverse) of XWX and multiplying on
  // the right by XWZ, i.e. B = (XWX)^-1 * XWZ. This new set of coefficients is
  // then compared to the current set (beta) and checked for convergence.
  bool ShouldIterate(<?=$constantState?>& modible_state) {
<?  if ($isGuassian) { ?>
    modible_state.iteration++;
    modible_state.beta = pinv(XWX, <?=$epsilon?>) * XWZ;
    return modible_state.iteration < 2;
<?  } else {
        if ($debug) { ?>
    // Debug Statements
    cout << "iteration: " << modible_state.iteration << endl;
    cout << "Performing calculation" << endl;
    cout << "XWX:" << endl << XWX << endl;
    cout << "XWZ:" << endl << XWZ << endl;
<?      } ?>
    modible_state.iteration++;
    // This return is done manually to avoid using pinv on 0 matrices. The XWX
    // and XWZ matrices are 0 because they aren't updated in the final run.
    // Computation of the new coefficients.
    <? /*Old version: vec new_beta = pinv(XWX, <?=$epsilon?>) * XWZ;*/?>
    vec new_beta = solve(XWX, XWZ);
    // TODO: Check that R uses this branch
    bool has_converged = modible_state.iteration == 1
      ? false
      : HasConverged(modible_state.beta, new_beta);
    modible_state.beta = new_beta;

    // For the iteration to conclude, either the max iteration must have been
    // reached or the new set of coefficients is sufficiently close to the
    // previous iteration. The last check is to ensure fitting was actually
    // performed, i.e. this wasn't the iteration computing starting values.
    return    modible_state.iteration < <?=$maxIteration, PHP_EOL?>
           && !has_converged;
<?  } ?>
  }

  bool HasConverged(vec old_beta, vec new_beta) {
    return max(abs((old_beta - new_beta) / old_beta)) < <?=$epsilon?>;
  }

  void GetResult(<?=typed_ref_args($outputs)?>) {
<?= ProduceResult(array_keys($outputs)[0], [
	'iteration' => 'constant_state.iteration',
	'coefficients' => 'constant_state.beta',
	'deviance' => 'deviance',
	'pearson' => 'pearson' ]) ?>

<?  if ($debug) { ?>
	cout << <?=array_keys($outputs)[0]?>.get() << endl;
<? } ?>
  }

  vec::fixed<<?=$numberCoefficients?>> GetCoefficients() const {
    return constant_state.beta;
  }
};

<?  return array(
        'kind' => 'GLA',
        'name' => $className,
        'system_headers' => [
            'math.h',
            'armadillo',
            'iostream',
            'boost/math/constants/constants.hpp',
            'boost/math/special_functions/erf.hpp',
        ],
        'lib_headers' => ['glmFamilies', 'ArmaJson'],
        'user_headers' => ['json.h'],
        'iterable' => TRUE,
        'input' => $inputs,
        'output' => $outputs,
        'result_type' => 'single',
        'generated_state' => $constantState,
        'libraries' => ['armadillo'],
    );
}

// This function should be called immediately on the terms object of the model.
// It bridges the interpration of R to the PHP GLA and processes the information
// taken from R. This is because R can interpret the model terms but it does not
// process them and has this do most of the processing, allowing it to be more
// diverse and simplifying the process. Also, in the case of factors, the
// cardinatilities are not known to R as the R interface is separated from the
// data.
// Inputs is an associative array with variable names pointing to type.
// Predictors is the terms object of the model, i.e. the array of predictors.
// Order is an array that specifies the order of each element of predictors.
function Construct_Input_Terms(array $inputs, array $predictors, array $order)
{
    $terms = array();
    for ($predictorCounter = 0;
         $predictorCounter < count($predictors); $predictorCounter ++) {
        // A list of the interacting terms is made
        $interactions = explode(':', $predictors[$predictorCounter]);
        // In R, you can specify a formula that has interactions in its AsIs
        // functions, such as 'y ~ I(x + w:z)`, which isn't allowed. This does
        // not result in an error when the formula is parsed but rather when
        // the model is constructed. Hence, this GLA needs to check for the
        // error itself rather than relying checks elsewhere.
        grokit_assert(
            count($interactions) == $order[$predictorCounter],
            'Improper formula - interactions were in an AsIs function.'
        );
        // Interacting terms that are numeric. Only contains names, as the
        // discrepancy between types (e.g. BIGINT vs DOUBLE) is unimportant.
        // All numeric types are treated as DOUBLE.
        $numericTerms  = array();
        // Interacting terms that are a factor. Names point to cardinality.
        $factorTerms   = array();
        for ($interactionCounter = 0;
             $interactionCounter < count($interactions);
             $interactionCounter ++) {
            $currentInteraction = $interactions[$interactionCounter];
            // Strip outside parenthesis from the interaction, most likely due
            // to AsIs functions. R allows you to specify a term with an
            // arbitrary arrangement of parenthesis. For example, 'y ~ (((x)))'
            // is valid. This would interfere in the type look-up below. Hence
            // this is needed both to clean up the output and to safeguard
            // against poorly specified predictor names.
            // Edit: This is properly redundant due to the use of 'format' in R.
            while ($currentInteraction[0] == "(")
                if (substr($currentInteraction, -1) == ")")
                    $currentInteraction =
                        substr($currentInteraction, 1, -1);
                else
                    break;
            // The variable is looked up in $inputs to learn its type..
            if (array_key_exists($currentInteraction, $inputs))
                if ($inputs[$currentInteraction]->is("numeric"))
                    $numericTerms[] = $currentInteraction;
                else
                    $factorTerms[$currentInteraction]
                        = $inputs[$currentInteraction]->get("cardinality");
            // The look-up failed and this interaction term is an AsIs function.
            // All AsIs functions are numeric.
            // TODO: Make sure there are no factors in the AsIs function.
            else
                $numericTerms[] = $currentInteraction;
        }
        // The arrays are sorted to ensure uniformity of output and to help with
        // the function Remove_Redundant_Keys
        sort($numericTerms);
        ksort($factorTerms);
        $terms[] = array($numericTerms, $factorTerms);
    }
    return $terms;
}

// This function should be called on the output of constructInputTerms. It
// removes any terms that are redundant due to the interaction of factors. For
// example, 'a' is redundant in 'a:b + a' where both 'a' and 'b' are factors. It
// does not remove terms that are redundant due to an AnIs statement in the
// model. For example, it would not remove anything from 'I(x*y) + x:y' where
// 'x' and 'y' are numeric. It operates on the reference of an array and hence
// does not return anything. The sorting performed in Construct_Input_Terms
// should be noted. This ensures the equality operator is viable.
// TODO: Be more explicit about what individual parts of this do and the structure of input.
function Remove_Redundant_Keys(array &$terms)
{
    // The counter iterates top down because higher order terms appear last.
    for ($termCounter= count($terms) - 1; $termCounter > 0; $termCounter --)
        for ($nestedCounter = $termCounter - 1; $nestedCounter >= 0;
             $nestedCounter --) {
            // As mentioned above, the elements of $terms are sorted.
            // Hence the array equality in the third condition is valid.
            $dummyArray1 = array_keys($terms[$nestedCounter][1]);
            $dummyArray2 = array_keys($terms[$termCounter][1]);
            if (   $termCounter != $nestedCounter
                && count(array_diff($dummyArray1, $dummyArray2)) == 0
                && $terms[$nestedCounter][0] == $terms[$termCounter][0]
            ) {
                // This ensures the removal of terms does not outpace the
                // decrement of termCounter so an out of bounds exception
                // doesn't occur.
                $termCounter --;
                array_splice($terms, $nestedCounter, 1);
            }
        }
}

// This is generates a power set of a set. It was taken from stack overflow and
// is likely to be correct. It is not cleaned up because it is most likely will
// not ever be called.
// TODO: Either delete this or implement it and clean up the code.

function Power_Set(array $in,$minLength = 1)
{
    $count = count($in);
    $members = pow(2,$count);
    $return = array();
    for ($i = 0; $i < $members; $i++) {
        $b = sprintf("%0".$count."b",$i);
        $out = array();
        for ($j = 0; $j < $count; $j++) {
            if ($b{$j} == '1') $out[] = $in[$j];
        }
        if (count($out) >= $minLength) {
          $return[] = $out;
        }
    }
    return $return;
}

// This generates a power set of an input set using recursion. Originally, it
// was needed for ANOVA. It is never called because the ANOVA is not implemented
// due to inefficiency.
// TODO: Implement this (and ANOVA in general) and test it.
function Power_Set_Recursive(array $inputSet)
{
    if(count($inputSet) == 1)
        return array($inputSet);
    $lastElement = array_pop($inputSet);
    $recursivePowerSet = powerSetRecursive($inputSet);
    $duplicateFamily = array();
    foreach ($recursivePowerSet as $currentSet)
        array_push($duplicateFamily,
                   array_merge($currentSet, array($lastElement)));
    return array_merge($recursivePowerSet, $duplicateFamily);
}

// This generates a family of sequential sets. Each set is the previous set plus
// another element in the input set, in order. For example, {a, b, c} would
// produce {{}, {a}, {a, b}, {a, b, c}}
// TODO: Implement this with ANOVA and test it.
function Sequential_Set(array $inputSet)
{
    $outputFamily = array(array());
    $currentSet = array();
    while(!empty($inputSet)) {
        $currentSet[]    =  array_shift($inputSet);
        $outputFamily[] =  $currentSet;
    }
    return $outputSet;
}

// Outputs c++ code to fill the vector that is constructed by an item. The
// vector's index for this term is equal to a function of its factors + the
// offset. The offset is considered in writeInputVector and represents the part
// of the vector already filled by previous terms. It is 0 if there are no
// factors. The exact function that specifies the index is that the index is the
// sum of the product of the value of each factor and the cardinalities of the
// factors after it. For example, in 'a:b:c' with |a| = 4, |b| = 3, |c| = 2, the
// index is 'a * 6 + b * 2 + c'. It should be noted that values for factors
// begin at 0 and hence the above expression ranges from 0 to 23, which agrees
// with c++ arrays beginning at index 0.
function Factors_To_Index(array $factors) {
    if(empty($factors))
        return "";
    $output = '';
    for($factorCounter = 0; $factorCounter < count($factors); $factorCounter ++)
        $output .= ' + ' . array_keys($factors)[$factorCounter] . '.GetID() * '
                .  array_product(array_slice($factors, $factorCounter + 1));
    return $output;
}

// Outputs c++ code to fill the vector that is constructed by an item. The
// vector is filled with the product of all numeric interactions for that term.
// Otherwise, it is set to 1 because that term only contains interactions
// between factor(s).
function Numeric_To_Value(array $numericInteractions) {
    if(empty($numericInteractions))
        return '1';
    else
        // Parenthesis included to ensure correct order of operations for AsIs.
        return '(' . implode(') * (', $numericInteractions) . ')';
}

// The following creates the code to construct an armadillo vector for the data
// item. The data needs to be packaged into a vector to utilize armadillo's
// speed. The return is array containing the size of the vector and the code to
// fill it with the item's data. Level is used to specify the level of
// identation and should usually be 3.
// TODO: Ensure level is correct.
function writeInputVector(array $terms, $interceptNeeded) {
    $offset = (int) $interceptNeeded;  // the offset of the c++ vector
    $cppCode = array();  // lines of code to be printed in the c++ file
    if ($interceptNeeded)
        $cppCode[] = 'x[0] = 1;';
    for ($termCounter = 0; $termCounter < count($terms); $termCounter ++) {
        $cppCode[] = "x[$offset"
                   . Factors_To_Index($terms[$termCounter][1]) . '] = '
                   . Numeric_To_Value($terms[$termCounter][0]) . ';';
        $offset += array_product($terms[$termCounter][1]);
    }
    // The inner implode creates the number of tabs specified by level.
    return array($offset, implode(PHP_EOL, $cppCode));
}

?>
