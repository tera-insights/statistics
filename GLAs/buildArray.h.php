<?
// This GLA is currently unfished.

// This GLA attempts to abstract the process of incrementally constructing multi
// dimensional arrays as much as possible, including the format of the input and
// the dimension of the array being made.
//
// The arguments should be a list of scalars followed by an array. The result is
// an array whose dimension is the number of scalars plus the dimension of the
// input array. Essentially, the list of scalars specify the coordinates of the
// input array in the final result. Each input should be numeric or categorical.
// The inputs specify the slice, column, and row index in that order. An example
// is that the input (integer, vector) specifies a matrix, in which the integer
// is the column index of the corresponding vector. Likewise, (integer, matrix)
// results in a cube and the integer gives the slice index of the matrix.
//
// Currently, only the construction of up to three-dimensional arrays (cubes) is
// supported, as that is what armadillo currently supports. This is unlikely to
// change except along with updates to armadillo to include higher dimensions of
// arrays.
//
// Template arguments:
// Aggregate: How to aggregate data when a single point in the result is given
//   by multiple indices. The default is sum, with multiplication, replacement
//   and omission also being supported. The use of replacement and omission is
//   not recommended, as they are not commutative and depend on the order of the
//   input.
//
// Allocation: An array with one element per dimension of the array being built.
//   Each element should be either a scalar or an array of scalars. In the first
//   case, the single scalar is simply converted to an array containing only it.
//   The first element of each such array should be a string specifying the
//   mode of input and allocation of memory for the result. Additional elements
//   parametrize that strategy. Each parameter should be named in the array and
//   the mode should never be named. The various modes are listed below:
//
//   Default: Typically the most efficient and most sensible solution. If the
//     corresponding input is a categorical variable, its cardinality is used as
//     the size for the sparse mode. Otherwise the sparse mode is used.
//     Additional parameters are the same for the designated mode, with a
//     warning being thrown if parameters are given for the other mode. If the
//     size parameter is given, it overrides the cardinality being used as such.
//
//   Sparse: A sparse array is built by storing both the indices and the data
//     separately. Indices marked as sparse are stored together separately from
//     the array built using non-sparse indices. See the Sparse parameter for
//     more information.
//
//   Fixed: The size of the result is fixed before any inputs are given. The
//     size parameter should be given as a positive integer. The check parameter
//     should be a string that matches "none", "error", or "ignore". If check is
//     "none", then no check is done to ensure the input is valid and in the
//     interval of [0, size). As such, this strategy is quite risky and at best,
//     a segmentation is thrown when bad input is given. If check is "error", an
//     error is thrown when bad input is first encountered. Lastly, when check
//     is "ignore", bad input is simply ignored. By default, check is "none" as
//     that is the most efficient option.
//
//   Sufficient: The structure is enlarged to accomodate new input. Additional
//     parameters are initial and maximum, which are both expected to either be
//     null or a positive integer. Initial specifies the initial size of the
//     allocated structure, with null denoting that it is initially empty. The
//     maximum parameter specifies the maximum allowed size. If maximum is not
//     null, the parameter error is allowed, which should be a boolean. If true,
//     an error is thrown when the input exceeds the maximum allowed. Otherwise,
//     it is simply ignored. Note that if that an index whose value is exactly
//     maximum is still barred, as indexing starts at 0. This strategy should be
//     used with some care as it can result memory overflow. However, it can be
//     quite useful with dealing with inputs of unknown size. Simply ordering
//     the input beforehand in descending order ensures that no-reallocation is
//     performed and the result is exactly the necessary size. It is also highly
//     recommended that only the outermost non-sparse dimension uses this mode.
//
//   Due the layout of memory used by the structure, the allocation modes cannot
//   be used in any order. Fixed and Sufficient cannot be used before Sparse due
//   to memory layout. Default attempts to work around with this. One example is
//   that (Default, Sparse, Fixed) is allowed even when the first input is a
//   categorical variable. While normally this would cause the first mode to be
//   Fixed, it is instead treated as Sparse. Additionally, there are several
//   workarounds to these restrictions. For example, (Fixed, Sparse) behaviour
//   can be replicated by using a GroupBy on the Fixed input with the inner GLA
//   being simply Sparse.
//
//   In addition to be the aforementioned array, Allocation can be a string that
//   specifies a mode, in which case that mode is used for each input using the
//   default parameters for that mode. If Allocation is given as an array, each
//   scalar argument can be a positive integer instead of a string, meaning that
//   the Fixed mode is used for that input with size equal to that integer.
//
// Sparse: An associative array of parameters governing the inner workings of
//   the sparse matrix allocation. A warning is given if this is specified when
//   no sparse modes are given.
//
//   Hash: A boolean denoting whether the inputs should be stored using a hash
//     mapping or a simple vector. In the former case, a hash is used to map a
//     tuple of sparse indices to the location of the corresponding data in the
//     structure used to hold non-sparse data. If false, the tuples are simply
//     kept in a vector, making it much harder to perform indexing on the data.
//     Furthermore, elements of repeated indices are aggregated only if hash is
//     true. Otherwise, repeated indices simply result in duplicated values in
//     the vector of tuples. The default value is false as is the more efficient
//     option when it comes to construction.
//
//   Size: The initial size of the object used to hold the tuples. This should
//     be a positive integer. The default value is NULL, meaning the default of
//     the associated data structure is used.
//
// Lastly, it should be noted that sparse matrices circumvent the limit of the
// dimension of the result. Although the result can only have dimension of at
// most 3, a sparse matrix is not subject to this as only the Fixed inputs must
// obey the limit. It is acceptable to have (Sparse, Sparse, Fixed) followed by
// a matrix input even though this simulates a result of dimension 5. However,
// input such as (Sparse, Fixed, Fixed) followed by a matrix input is still not
// allowed, as the Fixed inputs plus the dimension of the last input is 4.
function Build_Array(array $t_args, array $inputs, array $outputs)
{
    // Class name is randomly generated.
    $className = generate_name("BuildArray");

    // The number of scalar indices provided.
    $numScalars = count($inputs) - 1;

    // At least two arguments must be provided, as using only one argument would
    // be no better than a simple aggregate on the single input.
    grokit_assert($numScalars >= 1,
                  'Build Array: at least 2 inputs expected.');

    $indices = array_slice($inputs, 0, $numScalars, true);
    // Each index is expected to be a scalar. Currently, no scalar property is
    // implemented so each input is simply ensured to be categorical or numeric.
    foreach ($indices as $name => $type)
        grokit_assert($type->is('numeric') || $type->is('categorical'),
                      "Build Array: $name has illegal type $type.");
    $indices = array_values($indices);

    // The aggregate is grabbed and checked for correctness.
    $aggregates = ['sum', 'multiplication', 'replacement', 'omission'];
    $aggregate = strtolower(get_default($t_args, 'aggregate', 'sum'));
    grokit_assert(in_array($aggregate, $aggregates),
                  "Build Array: illegal aggregate: $aggregate");

    // The allocation strategy is grabbed and formatted correctly.
    $allocation = array_values(get_default($t_args, 'allocation', 'default'));
    if (is_string($allocation))
        $allocation = array_fill(0, $numScalars, $allocation);
    $type = gettype($allocation);
    grokit_assert($type == 'array',
                  "Build Matrix: illegal type for allocation: $type");
    grokit_assert(count($allocation) == $numScalars,
                  'Build Matrix: incorrect number of elements in allocation.');
    foreach ($allocation as &$element) {
        $type = gettype($element)
        grokit_assert(in_array($type, ['array', 'string', 'integer']),
                      "Build Matrix: illegal type for allocation entry: $type");
        if ($type == 'string')
            $element = [$element];
        else if ($type == 'integer')
            $element = ['fixed', 'size' => $element];
        else if ($type == 'array')
            $element = array_change_key_case($element, CASE_LOWER);
        else
            grokit_error("Build Matrix: illegal allocation entry type: $type");
    }

    // Separating the modes and parameters for each allocation strategy.
    foreach ($allocation as $element) {
        $keys = array_keys($element);
        $good = !in_array(false, array_map('is_string', array_slice($keys, 1)))
              & is_integer($keys[0])
              & $keys[0] == 0;
        grokit_assert($good, 'Build Matrix: invalid keys in allocation entry.');
        $modes[] = strtolower($element[0]);
        $param[] = array_slice($element, 1, NULL, true);
    }

    $checks = [
        'sparse' => [],
        'fixed' => [
            'size'  => ['is_integer', 'is_positive'],
            'error' => ['is_string', 'in_array', [['none', 'error', 'ignore']]],
        ],
        'sufficient' => [
            'initial' => ['or_func', [[['is_null'], ['is_integer', 'is_positive']]]],
            'maximum' => ['or_func', [[['is_null'], ['is_integer', 'is_positive']]]],
            'error'   => ['is_bool'];
        ],
    ];

    $defaults = [
        'sparse' => [],
        'fixed' => [
            'error' => 'none',
        ],
        'sufficient' => [
            'initial' => NULL,
            'maximum' => NULL,
            'error'   => false,
        ],
    ];

    // Validation of modes.
    $valid = ['default', 'sparse', 'fixed', 'sufficient'];
    $invalid = array_diff($modes, $valid);
    grokit_assert(count($invalid) == 0,
                  'Build Array: invalid modes: ' . implode(', ', $invalid));

    // Decoding of default types into sparse and fixed is done here.
    $key = in_array('sparse', $modes) ? array_search('sparse', $modes) : -1;

    for ($i = 0; $i < $key; $i++)
        if ($modes[$i] == 'default')
            $modes[$i] = 'sparse';
        else
            grokit_assert($modes[$i] == 'sparse',
                          'Build Array: mode $modes[$i] used before sparse.');

    for ($i = $key + 1; $i < $numScalars; $i++) {
        if ($modes[$i] == 'default') {
            if ($indices[$i]->is('categorical')) {
                $modes[$i] = 'fixed';
                if (!array_key_exists('size', $param[$i]))
                    $param[$i]['size'] = $indices[$i]->get('cardinality');
            } else {
                $modes[$i] = 'sparse';
            }
        }
    }

    // Each allocation mode now has its parameters filled in and verified.
    $strategies = array_combine($modes, $param);
    foreach ($strategies as $mode => &$args)
        process_args($args, $defaults[$mode], $checks[$mode]);

    $isSparse = in_array('sparse', $modes);
    $sparse = get_default($t_args, 'sparse', NULL);
    grokit_warning_if(!($isSparse || is_null($sparse)),
                      'Build Array: sparse parameters given unnecessarily.');

    $check = [
        'size' => ['or_func', [[['is_null'], ['is_integer', 'is_positive']]]],
        'hash' => ['is_bool'],
    ];
    $default = [
        'size' => null,
        'hash' => false,
    ];

    if (is_null($sparse))
        $sparse = $default;

    grokit_assert(is_array($sparse),
                  'Build Array: sparse should be null or an array.');

    process_args($sparse, $default, $check);

    // Whether all dimensions are fixed size.
    $isFixed = array_unique($modes) == ['fixed'];

    // The type of the data being assimiliated.
    $type = array_get_index($inputs, $numScalars);

    // The dimension of the result.
    $dimension = $numScalars - $key - !$isSparse + get_dimension($type);
}
?>