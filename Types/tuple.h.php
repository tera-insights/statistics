<?
//  Copyright 2014 Tera Insights, LLC. All Rights Reserved.

// An stl tuple container.

function Tuple(array $t_args) {
    $types = $t_args['types'];

    grokit_assert(is_array($types), 'Tuple: [types] should be an array.');
    foreach ($types as $type)
        grokit_assert(is_datatype($value), "Tuple: [$type] is not a datatype.");

    $className = generate_name('Tuple');

    $sys_headers     = ['tuple'];
    $user_headers    = [];
    $lib_headers     = [];
    $constructors    = [];
    $methods         = [];
    $functions       = [];
    $binaryOperators = [];
    $unaryOperators  = [];
    $globalContent   = '';
    $complex         = false;
    $properties      = ['tuple'];
    $extras          = ['types' => $types, 'size' => count($types)];
?>

using std::tuple<<?=typed($types)?>> = <?=$className?>;

<?
    return [
        'kind'             => 'TYPE',
        'name'             => $className,
        'system_headers'   => $sys_headers,
        'user_headers'     => $user_headers,
        'lib_headers'      => $lib_headers,
        'constructors'     => $constructors,
        'methods'          => $methods,
        'functions'        => $functions,
        'binary_operators' => $binaryOperators,
        'unary_operators'  => $unaryOperators,
        'global_content'   => $globalContent,
        'complex'          => $complex,
        'properties'       => $properties,
        'extras'           => $extras,
    ];
}
?>
