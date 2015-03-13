<?
//  Copyright 2014 Tera Insights, LLC. All Rights Reserved.

// A fixed size matrix from Armadillo.

function Fixed_Vector($t_args) {
    $size = get_default($t_args, 'size', null);
    $direction = get_default($t_args, 'direction', 'col');
    $type = lookupType(get_default($t_args, 'type', "base::double"));

    grokit_assert(is_string($direction) && in_array($direction, ['row', 'col']),
                  'Vector: [direction] argument must be "row" or "col".');
    grokit_assert(is_datatype($type) && $type->is('numeric'),
                  'Vector: [type] argument must be a numeric datatype.');
    grokit_assert(is_int($size) && $size > 0,
                  'Vector: [size] must be a positive integer.');

    $inputs = get_default($t_args, 'inputs', array_fill(0, $size, $type));
    grokit_assert(is_array($inputs),
                  'Vector: [inputs] must be an array of datatypes.');

    foreach ($inputs as $key => &$input) {
        grokit_assert(   is_datatype($input)
                      && ($input->is('numeric') || $input->is('categorical')),
                      "MakeVector: input [$key] must be a numeric type.");
        $input = $input->lookup();
    };

    $className = generate_name('Vector_' . $size . '_');

    $cppType = ($direction == 'row') ? 'Row' : 'Col';

    $unique = ['Vector', $direction, $type, $size];
    if (typeDefined($unique, $hash, $output))
      return $output;

    $sys_headers     = ['armadillo'];
    $user_headers    = [];
    $lib_headers     = [];
    $constructors    = [];
    $methods         = [];
    $functions       = [];
    $binaryOperators = [];
    $unaryOperators  = [];
    $globalContent   = '';
    $complex         = false;
    $properties      = ['vector'];
    $extras          = ['size' => $size, 'direction' => $direction,
                        'type' => $type, 'inputs' => $inputs];
?>

typedef arma::<?=$cppType?><<?=$type?>>::fixed<<?=$size?>> <?=$className?>;

<?  ob_start(); ?>

inline void FromString(@type& x, const char* buffer) {
  char * current = NULL;
  char * saveptr = NULL;
  const char * delim = " ";
  char * copy = strdup(buffer);

  current = strtok_r(copy, delim, &saveptr);

  for (auto& val : x) {
    FATALIF(current == NULL, "Not enough elements in string representation of array");
    ToString(val, current);
    current = strtok_r(NULL, delim, &saveptr);
  }

  free((void *) copy);
}

inline void ToJson(const @type src, Json::Value& dest) {
  dest["__type__"] = "vector";
  dest["n_elem"] = src.n_elem;
  Json::Value content(Json::arrayValue);
  for (int i = 0; i < src.n_elem; i++)
    content[i] = src(i);
  dest["data"] = content;
}

inline size_t Serialize(char* buffer, const @type& src) {
  <?=$type?> * asInnerType = reinterpret_cast<<?=$type?>*>(buffer);
  const <?=$type?> * colPtr = src.memptr();
  std::copy(colPtr, colPtr + @type::n_elem, asInnerType);
  return @type::n_elem * sizeof(<?=$type?>);
}

inline size_t Deserialize(const char* buffer, @type& dest) {
  const <?=$type?>* asInnerType = reinterpret_cast<const <?=$type?>*>(buffer);
  dest = @type(asInnerType);
  return @type::n_elem * sizeof(<?=$type?>);
}

inline size_t SerializedSize(@type& dest) {
  return @type::n_elem * sizeof(<?=$type?>);
}

<?  $globalContent .= ob_get_clean(); ?>

<?
    return [
        'kind'             => 'TYPE',
        'name'             => $className,
        'hash'             => $hash,
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

declareType('Vector', 'statistics::Fixed_Vector', []);
declareType('FixedVector', 'statistics::Fixed_Vector', []);
?>
