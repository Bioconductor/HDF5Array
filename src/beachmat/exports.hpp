#define REGISTER(x) R_RegisterCCallable("HDF5Array", #x, reinterpret_cast<DL_FUNC>(x))

REGISTER(HDF5Matrix_integer_input_create);

REGISTER(HDF5Matrix_integer_input_destroy);

REGISTER(HDF5Matrix_integer_input_clone);

REGISTER(HDF5Matrix_integer_input_dim);

REGISTER(HDF5Matrix_integer_input_get);

REGISTER(HDF5Matrix_integer_input_getRow_integer);

REGISTER(HDF5Matrix_integer_input_getCol_integer);

REGISTER(HDF5Matrix_integer_input_getRow_numeric);

REGISTER(HDF5Matrix_integer_input_getCol_numeric);

REGISTER(HDF5Matrix_integer_input_getRows_integer);

REGISTER(HDF5Matrix_integer_input_getCols_integer);

REGISTER(HDF5Matrix_integer_input_getRows_numeric);

REGISTER(HDF5Matrix_integer_input_getCols_numeric);

REGISTER(HDF5Matrix_integer_output_create);

REGISTER(HDF5Matrix_integer_output_destroy);

REGISTER(HDF5Matrix_integer_output_clone);

REGISTER(HDF5Matrix_integer_output_yield);

REGISTER(HDF5Matrix_integer_output_get);

REGISTER(HDF5Matrix_integer_output_getRow_integer);

REGISTER(HDF5Matrix_integer_output_getCol_integer);

REGISTER(HDF5Matrix_integer_output_getRow_numeric);

REGISTER(HDF5Matrix_integer_output_getCol_numeric);

REGISTER(HDF5Matrix_integer_output_set);

REGISTER(HDF5Matrix_integer_output_setRow_integer);

REGISTER(HDF5Matrix_integer_output_setCol_integer);

REGISTER(HDF5Matrix_integer_output_setRow_numeric);

REGISTER(HDF5Matrix_integer_output_setCol_numeric);

REGISTER(HDF5Matrix_integer_output_setRowIndexed_integer);

REGISTER(HDF5Matrix_integer_output_setColIndexed_integer);

REGISTER(HDF5Matrix_integer_output_setRowIndexed_numeric);

REGISTER(HDF5Matrix_integer_output_setColIndexed_numeric);

REGISTER(HDF5Matrix_logical_input_create);

REGISTER(HDF5Matrix_logical_input_destroy);

REGISTER(HDF5Matrix_logical_input_clone);

REGISTER(HDF5Matrix_logical_input_dim);

REGISTER(HDF5Matrix_logical_input_get);

REGISTER(HDF5Matrix_logical_input_getRow_integer);

REGISTER(HDF5Matrix_logical_input_getCol_integer);

REGISTER(HDF5Matrix_logical_input_getRow_numeric);

REGISTER(HDF5Matrix_logical_input_getCol_numeric);

REGISTER(HDF5Matrix_logical_input_getRows_integer);

REGISTER(HDF5Matrix_logical_input_getCols_integer);

REGISTER(HDF5Matrix_logical_input_getRows_numeric);

REGISTER(HDF5Matrix_logical_input_getCols_numeric);

REGISTER(HDF5Matrix_logical_output_create);

REGISTER(HDF5Matrix_logical_output_destroy);

REGISTER(HDF5Matrix_logical_output_clone);

REGISTER(HDF5Matrix_logical_output_yield);

REGISTER(HDF5Matrix_logical_output_get);

REGISTER(HDF5Matrix_logical_output_getRow_integer);

REGISTER(HDF5Matrix_logical_output_getCol_integer);

REGISTER(HDF5Matrix_logical_output_getRow_numeric);

REGISTER(HDF5Matrix_logical_output_getCol_numeric);

REGISTER(HDF5Matrix_logical_output_set);

REGISTER(HDF5Matrix_logical_output_setRow_integer);

REGISTER(HDF5Matrix_logical_output_setCol_integer);

REGISTER(HDF5Matrix_logical_output_setRow_numeric);

REGISTER(HDF5Matrix_logical_output_setCol_numeric);

REGISTER(HDF5Matrix_logical_output_setRowIndexed_integer);

REGISTER(HDF5Matrix_logical_output_setColIndexed_integer);

REGISTER(HDF5Matrix_logical_output_setRowIndexed_numeric);

REGISTER(HDF5Matrix_logical_output_setColIndexed_numeric);

REGISTER(HDF5Matrix_numeric_input_create);

REGISTER(HDF5Matrix_numeric_input_destroy);

REGISTER(HDF5Matrix_numeric_input_clone);

REGISTER(HDF5Matrix_numeric_input_dim);

REGISTER(HDF5Matrix_numeric_input_get);

REGISTER(HDF5Matrix_numeric_input_getRow_integer);

REGISTER(HDF5Matrix_numeric_input_getCol_integer);

REGISTER(HDF5Matrix_numeric_input_getRow_numeric);

REGISTER(HDF5Matrix_numeric_input_getCol_numeric);

REGISTER(HDF5Matrix_numeric_input_getRows_integer);

REGISTER(HDF5Matrix_numeric_input_getCols_integer);

REGISTER(HDF5Matrix_numeric_input_getRows_numeric);

REGISTER(HDF5Matrix_numeric_input_getCols_numeric);

REGISTER(HDF5Matrix_numeric_output_create);

REGISTER(HDF5Matrix_numeric_output_destroy);

REGISTER(HDF5Matrix_numeric_output_clone);

REGISTER(HDF5Matrix_numeric_output_yield);

REGISTER(HDF5Matrix_numeric_output_get);

REGISTER(HDF5Matrix_numeric_output_getRow_integer);

REGISTER(HDF5Matrix_numeric_output_getCol_integer);

REGISTER(HDF5Matrix_numeric_output_getRow_numeric);

REGISTER(HDF5Matrix_numeric_output_getCol_numeric);

REGISTER(HDF5Matrix_numeric_output_set);

REGISTER(HDF5Matrix_numeric_output_setRow_integer);

REGISTER(HDF5Matrix_numeric_output_setCol_integer);

REGISTER(HDF5Matrix_numeric_output_setRow_numeric);

REGISTER(HDF5Matrix_numeric_output_setCol_numeric);

REGISTER(HDF5Matrix_numeric_output_setRowIndexed_integer);

REGISTER(HDF5Matrix_numeric_output_setColIndexed_integer);

REGISTER(HDF5Matrix_numeric_output_setRowIndexed_numeric);

REGISTER(HDF5Matrix_numeric_output_setColIndexed_numeric);
