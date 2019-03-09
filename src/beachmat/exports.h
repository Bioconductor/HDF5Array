#include "Rcpp.h"

extern "C" {

void * HDF5Matrix_integer_input_create (SEXP);

void HDF5Matrix_integer_input_destroy (void *);

void * HDF5Matrix_integer_input_clone (void *);

void HDF5Matrix_integer_input_dim(void*, size_t*, size_t*);
 
void HDF5Matrix_integer_input_get(void *, size_t, size_t, int*);

void HDF5Matrix_integer_input_getRow_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_input_getCol_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_input_getRow_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_input_getCol_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_input_getRows_integer(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_input_getCols_integer(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_input_getRows_numeric(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_input_getCols_numeric(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void * HDF5Matrix_integer_output_create (size_t, size_t);

void HDF5Matrix_integer_output_destroy (void *);

void * HDF5Matrix_integer_output_clone (void *);

SEXP HDF5Matrix_integer_output_yield (void *);

void HDF5Matrix_integer_output_get(void *, size_t, size_t, int*);

void HDF5Matrix_integer_output_getRow_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_output_getCol_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_output_getRow_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_output_getCol_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_output_set(void *, size_t, size_t, int*);

void HDF5Matrix_integer_output_setRow_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_output_setCol_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_output_setRow_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_output_setCol_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void HDF5Matrix_integer_output_setRowIndexed_integer(void *, size_t, size_t, Rcpp::IntegerVector::iterator*, Rcpp::IntegerVector::iterator*);

void HDF5Matrix_integer_output_setColIndexed_integer(void *, size_t, size_t, Rcpp::IntegerVector::iterator*, Rcpp::IntegerVector::iterator*);

void HDF5Matrix_integer_output_setRowIndexed_numeric(void *, size_t, size_t, Rcpp::IntegerVector::iterator*, Rcpp::NumericVector::iterator*);

void HDF5Matrix_integer_output_setColIndexed_numeric(void *, size_t, size_t, Rcpp::IntegerVector::iterator*, Rcpp::NumericVector::iterator*);

}
