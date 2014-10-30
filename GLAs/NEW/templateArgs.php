<?
require "/opt/datapath/src/PHP/php/grokit_base.php";

function templateArgs($pattern, $glue)
{
    if (func_num_args() < 3)
        grokit_error("No arrays given");
    $numArrays = func_num_args() - 2;
    // Eval evaluates at local scope, can't use func_get_arg
    for ($counter = 1; $counter <= $numArrays; $counter++) {
        $dummyVar = func_get_arg($counter + 1);
        eval("\$array$counter = \$dummyVar;");
    }
    for ($counter = 1; $counter <= $numArrays; $counter ++)
        eval("if (!is_array(\$array$counter))
                  grokit_error(\"Non-array template args\");");
    $length = count($array1);
    for ($counter = 2; $counter <= $numArrays; $counter++) {
        eval("\$currentLength = count(\$array$counter);");
        if ($length != $currentLength)
            grokit_error("Unequal lengths of array arguments.");
    }
    $codeLines = array(); // array of c++ code
    for ($counter = 0; $counter < $length; $counter ++) {
        for ($setter= 1; $setter <= $numArrays; $setter++) {
            eval("\$key$setter = array_keys(\$array$setter)[$setter - 1];");
            eval("\$val$setter = \$array{$setter}[\$key{$setter}];");
        }
        eval("\$codeLines[] = $pattern;");
    }
    return implode($glue, $codeLines);
}
