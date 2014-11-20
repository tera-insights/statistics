<?
//  Copyright 2014 Tera Insights, LLC. All Rights Reserved.

// A fixed size matrix from Armadillo.

function Fixed_Matrix(array $t_args) {
    $nrow = $t_args['nrow'];
    $ncol = $t_args['ncol'];
    $type = get_default($t_args, 'type', lookupType("base::double"));

    grokit_assert(is_datatype($type),
                  'Matrix: [type] argument must be a valid datatype.');
    grokit_assert($type->is('numeric'),
                  'Matrix: [type] argument must be a numeric datatype.');
    grokit_assert(is_int($ncol) && $ncol > 0,
                  'Matrix: [ncol] must be a positiver integer.');
    grokit_assert(is_int($nrow) && $nrow > 0,
                  'Matrix: [nrow] must be a positive integer.');

    $className = generate_name('Matrix_' . $nrow . '_' . $ncol . '_');
    $cppType = $type->is('real') ? 'mat' : 'imat';

    $sys_headers = ['armadillo'];
    $user_headers = [];
    $lib_headers = ['ArmaJson'];
    $constructors = [];
    $methods = [];
    $functions = [];
    $binaryOperators = [];
    $unaryOperators = [];
    $globalContent = '';
    $complex = false;
    $properties = [];
    $extra = ['nrow' => $nrow, 'ncol' => $ncol, 'type' => $type];
?>

typedef arma::<?=$cppType?>::fixed<<?=$nrow?>, <?=$ncol?>> <?=$className?>;

<?  ob_start(); ?>

inline void ToJson(const @type src, Json::Value& dest) {
  dest["__type__"] = "matrix";
  dest["n_rows"] = src.n_rows;
  dest["n_cols"] = src.n_cols;
  Json::Value content(Json::arrayValue);
  for (int i = 0; i < src.n_rows; i++)
    for (int j = 0; j < src.n_cols; j++)
	    content[i * src.n_cols + j] = src(i, j);
  dest["data"] = content;
}

<?  $globalContent .= ob_get_clean(); ?>

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
        'extra'            => $extra,
    ];
}

declareType('Matrix', 'statistics::Fixed_Matrix', []);
declareType('FixedMatrix', 'statistics::Fixed_Vector', []);
?>
