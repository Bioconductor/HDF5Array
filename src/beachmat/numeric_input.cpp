#include "LIN_matrix.h"
#include "exports.h"

typedef HDF5_lin_matrix<double, REALSXP> HDF5NumMat;

// Constructor, destructors and clones.

void * HDF5Matrix_numeric_input_create (SEXP incoming) {
    return static_cast<void*>(new HDF5NumMat(incoming));
}

void HDF5Matrix_numeric_input_destroy (void * ptr) {
    delete static_cast<HDF5NumMat*>(ptr);
    return;
}

void * HDF5Matrix_numeric_input_clone (void * ptr) {
    HDF5NumMat* old=static_cast<HDF5NumMat*>(ptr);
    return static_cast<void*>(new HDF5NumMat(*old));
}

// Basic getters

void HDF5Matrix_numeric_input_dim(void* ptr, size_t* nr, size_t* nc){ 
    HDF5NumMat* thing=static_cast<HDF5NumMat*>(ptr);
    *nr=thing->get_nrow();
    *nc=thing->get_ncol();
    return;
}

void HDF5Matrix_numeric_input_get(void * ptr, size_t r, size_t c, double* val) {
    *val=static_cast<HDF5NumMat*>(ptr)->get(r, c);
    return;
}

void HDF5Matrix_numeric_input_getRow_integer(void * ptr, size_t r, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumMat*>(ptr)->get_row(r, *out, first, last);
    return;
}

void HDF5Matrix_numeric_input_getCol_integer(void * ptr, size_t c, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumMat*>(ptr)->get_col(c, *out, first, last);
    return;
}

void HDF5Matrix_numeric_input_getRow_numeric(void * ptr, size_t r, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumMat*>(ptr)->get_row(r, *out, first, last);
    return;
}

void HDF5Matrix_numeric_input_getCol_numeric(void * ptr, size_t c, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumMat*>(ptr)->get_col(c, *out, first, last);
    return;
}

// Multi getters

void HDF5Matrix_numeric_input_getRows_integer(void * ptr, Rcpp::IntegerVector::iterator* r, size_t n, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumMat*>(ptr)->get_rows(*r, n, *out, first, last);
    return;
}

void HDF5Matrix_numeric_input_getCols_integer(void * ptr, Rcpp::IntegerVector::iterator* c, size_t n, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumMat*>(ptr)->get_cols(*c, n, *out, first, last);
    return;
}

void HDF5Matrix_numeric_input_getRows_numeric(void * ptr, Rcpp::IntegerVector::iterator* r, size_t n, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumMat*>(ptr)->get_rows(*r, n, *out, first, last);
    return;
}

void HDF5Matrix_numeric_input_getCols_numeric(void * ptr, Rcpp::IntegerVector::iterator* c, size_t n, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumMat*>(ptr)->get_cols(*c, n, *out, first, last);
    return;
}
