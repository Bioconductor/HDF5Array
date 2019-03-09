#include "LIN_output.h"
#include "exports.h"

template<>
Rcpp::RObject HDF5_writer<double, REALSXP>::get_firstval() { 
    double first=0;
    if (get_nrow() && get_ncol()) {
        extract_one(0, 0, &first);
    }
    return Rcpp::NumericVector::create(first);
}

template<>
double HDF5_writer<double, REALSXP>::get_empty() { return 0; }

typedef HDF5_lin_output<double, Rcpp::NumericVector, REALSXP> HDF5NumOut;

// Constructor, destructors and clones.

void * HDF5Matrix_numeric_output_create (size_t nr, size_t nc) {
    return static_cast<void*>(new HDF5NumOut(nr, nc));
}

void HDF5Matrix_numeric_output_destroy (void * ptr) {
    delete static_cast<HDF5NumOut*>(ptr);
    return;
}

void * HDF5Matrix_numeric_output_clone (void * ptr) {
    HDF5NumOut* old=static_cast<HDF5NumOut*>(ptr);
    return static_cast<void*>(new HDF5NumOut(*old));
}

SEXP HDF5Matrix_numeric_output_yield (void * ptr) {
    return static_cast<HDF5NumOut*>(ptr)->yield();
}

// Basic getters

void HDF5Matrix_numeric_output_get(void * ptr, size_t r, size_t c, double* val) {
    *val=static_cast<HDF5NumOut*>(ptr)->get(r, c);	
    return;
}

void HDF5Matrix_numeric_output_getRow_integer(void * ptr, size_t r, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumOut*>(ptr)->get_row(r, *out, first, last);
    return;
}

void HDF5Matrix_numeric_output_getCol_integer(void * ptr, size_t c, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumOut*>(ptr)->get_col(c, *out, first, last);
    return;
}

void HDF5Matrix_numeric_output_getRow_numeric(void * ptr, size_t r, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumOut*>(ptr)->get_row(r, *out, first, last);
    return;
}

void HDF5Matrix_numeric_output_getCol_numeric(void * ptr, size_t c, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumOut*>(ptr)->get_col(c, *out, first, last);
    return;
}

// Basic setters

void HDF5Matrix_numeric_output_set(void * ptr, size_t r, size_t c, double* val) {
    static_cast<HDF5NumOut*>(ptr)->set(r, c, *val);	
    return;
}

void HDF5Matrix_numeric_output_setRow_integer(void * ptr, size_t r, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumOut*>(ptr)->set_row(r, *out, first, last);
    return;
}

void HDF5Matrix_numeric_output_setCol_integer(void * ptr, size_t c, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumOut*>(ptr)->set_col(c, *out, first, last);
    return;
}

void HDF5Matrix_numeric_output_setRow_numeric(void * ptr, size_t r, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumOut*>(ptr)->set_row(r, *out, first, last);
    return;
}

void HDF5Matrix_numeric_output_setCol_numeric(void * ptr, size_t c, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5NumOut*>(ptr)->set_col(c, *out, first, last);
    return;
}

// Indexed setters

void HDF5Matrix_numeric_output_setRowIndexed_integer(void * ptr, size_t r, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::IntegerVector::iterator* out) {
    static_cast<HDF5NumOut*>(ptr)->set_row_indexed(r, n, *idx, *out);
    return;
}

void HDF5Matrix_numeric_output_setColIndexed_integer(void * ptr, size_t c, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::IntegerVector::iterator* out) {
    static_cast<HDF5NumOut*>(ptr)->set_col_indexed(c, n, *idx, *out);
    return;
}

void HDF5Matrix_numeric_output_setRowIndexed_numeric(void * ptr, size_t r, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::NumericVector::iterator* out) {
    static_cast<HDF5NumOut*>(ptr)->set_row_indexed(r, n, *idx, *out);
    return;
}

void HDF5Matrix_numeric_output_setColIndexed_numeric(void * ptr, size_t c, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::NumericVector::iterator* out) {
    static_cast<HDF5NumOut*>(ptr)->set_col_indexed(c, n, *idx, *out);
    return;
}
