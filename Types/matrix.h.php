<?
//  Copyright 2014 Tera Insights, LLC. All Rights Reserved.

// A fixed size matrix from Armadillo.

function Fixed_Matrix(array $t_args) {
    $nrow = $t_args['nrow'];
    $ncol = $t_args['ncol'];
    $type = get_default($t_args, 'type', lookupType("base::double"));
    // TODO: Automatically look-up type.
    lookupType($type);

    grokit_assert(is_datatype($type),
                  'Matrix: [type] argument must be a valid datatype.');
    grokit_assert($type->is('numeric'),
                  'Matrix: [type] argument must be a numeric datatype.');
    grokit_assert(is_int($ncol) && $ncol > 0,
                  'Matrix: [ncol] must be a positiver integer.');
    grokit_assert(is_int($nrow) && $nrow > 0,
                  'Matrix: [nrow] must be a positive integer.');

    $nelem = $nrow * $ncol;

    $className = generate_name('Matrix_' . $nrow . '_' . $ncol . '_');

    $sys_headers     = ['armadillo', 'algorithm'];
    $user_headers    = [];
    $lib_headers     = ['ArmaJson'];
    $constructors    = [];
    $methods         = [];
    $functions       = [];
    $binaryOperators = [];
    $unaryOperators  = [];
    $globalContent   = '';
    $complex         = "ColumnIterator<@type, 0, sizeof({$type}) * {$nelem}>";
    $properties      = ['matrix', 'armadillo', 'fixed'];
    $extra           = ['size' => $nrow * $ncol, 'dimensions' => [$nrow, $ncol],
                        'nrow' => $nrow, 'ncol' => $ncol, 'type' => $type];
?>

typedef arma::Mat<<?=$type?>>::fixed<<?=$nrow?>, <?=$ncol?>> <?=$className?>;

<?  ob_start(); ?>

inline void ToJson(const @type src, Json::Value& dest) {
  FATALIF(src.n_elem > 4294967295L  // Maximum value for a uint32_t.
          "Error: attempt to serialize matrix with more than 2^31 - 1 elements")

  dest["__type__"] = "matrix";
  dest["n_rows"] = src.n_rows;
  dest["n_cols"] = src.n_cols;
  Json::Value content(Json::arrayValue);
  for (arma::uword i = 0; i != src.n_rows; i++) {
    for (arma::uword j = 0; j != src.n_cols; j++) {
      // Force index to be  Json::ArrayIndex to prevent the compiler from
      // getting confused at to which operator[] to use on content.
      Json::ArrayIndex position = (i * src.n_cols) + j;
      content[position] = src(i, j);
    }
  }
  dest["data"] = content;
}

template<>
inline size_t Serialize(char* buffer, const @type& src) {
  <?=$type?>* asInnerType = reinterpret_cast<<?=$type?>*>(buffer);
  const <?=$type?>* colPtr = src.memptr();
  std::copy(colPtr, colPtr + @type::n_elem, asInnerType);
  return @type::n_elem * sizeof(<?=$type?>);
}

template<>
inline size_t Deserialize(const char* buffer, @type& dest) {
  const <?=$type?>* asInnerType = reinterpret_cast<const <?=$type?>*>(buffer);
  dest = @type(asInnerType);
  return @type::n_elem * sizeof(<?=$type?>);
}

template<>
inline size_t SerializedSize(const @type& dest) {
  return @type::n_elem * sizeof(<?=$type?>);
}

<?  $globalContent .= ob_get_clean(); ?>

<?
    $innerDesc = function($var, $myType) use($type, $ncol, $nrow) {
        $describer = $type->describer('json');
?>
        <?=$var?>["n_cols"] = Json::Int64(<?=$ncol?>);
        <?=$var?>["n_rows"] = Json::Int64(<?=$nrow?>);
<?
        $innerVar = "{$var}[\"inner_type\"]";
        $describer($innerVar, $type);
    };

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
        'describe_json'    => DescribeJson('matrix', $innerDesc),
    ];
}

declareType('Matrix', 'statistics::Fixed_Matrix', []);
declareType('FixedMatrix', 'statistics::Fixed_Vector', []);
?>
