#include "LIN_output.h"
#include "exports.h"

template<>
Rcpp::RObject HDF5_writer<int, INTSXP>::get_firstval() { 
    int first=0;
    if (get_nrow() && get_ncol()) {
        extract_one(0, 0, &first);
    }
    return Rcpp::IntegerVector::create(first);
}

template<>
int HDF5_writer<int, INTSXP>::get_empty() { return 0; }

typedef HDF5_lin_output<int, Rcpp::IntegerVector, INTSXP> HDF5IntOut;

// Constructor, destructors and clones.

void * HDF5Matrix_integer_output_create (size_t nr, size_t nc) {
    return static_cast<void*>(new HDF5IntOut(nr, nc));
}

void HDF5Matrix_integer_output_destroy (void * ptr) {
    delete static_cast<HDF5IntOut*>(ptr);
    return;
}

void * HDF5Matrix_integer_output_clone (void * ptr) {
    HDF5IntOut* old=static_cast<HDF5IntOut*>(ptr);
    return static_cast<void*>(new HDF5IntOut(*old));
}

SEXP HDF5Matrix_integer_output_yield (void * ptr) {
    return static_cast<HDF5IntOut*>(ptr)->yield();
}

// Basic getters

void HDF5Matrix_integer_output_get(void * ptr, size_t r, size_t c, int* val) {
    *val=static_cast<HDF5IntOut*>(ptr)->get(r, c);	
    return;
}

void HDF5Matrix_integer_output_getRow_integer(void * ptr, size_t r, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5IntOut*>(ptr)->get_row(r, *out, first, last);
    return;
}

void HDF5Matrix_integer_output_getCol_integer(void * ptr, size_t c, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5IntOut*>(ptr)->get_col(c, *out, first, last);
    return;
}

void HDF5Matrix_integer_output_getRow_numeric(void * ptr, size_t r, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5IntOut*>(ptr)->get_row(r, *out, first, last);
    return;
}

void HDF5Matrix_integer_output_getCol_numeric(void * ptr, size_t c, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5IntOut*>(ptr)->get_col(c, *out, first, last);
    return;
}

// Basic setters

void HDF5Matrix_integer_output_set(void * ptr, size_t r, size_t c, int* val) {
    static_cast<HDF5IntOut*>(ptr)->set(r, c, *val);	
    return;
}

void HDF5Matrix_integer_output_setRow_integer(void * ptr, size_t r, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5IntOut*>(ptr)->set_row(r, *out, first, last);
    return;
}

void HDF5Matrix_integer_output_setCol_integer(void * ptr, size_t c, Rcpp::IntegerVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5IntOut*>(ptr)->set_col(c, *out, first, last);
    return;
}

void HDF5Matrix_integer_output_setRow_numeric(void * ptr, size_t r, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5IntOut*>(ptr)->set_row(r, *out, first, last);
    return;
}

void HDF5Matrix_integer_output_setCol_numeric(void * ptr, size_t c, Rcpp::NumericVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5IntOut*>(ptr)->set_col(c, *out, first, last);
    return;
}

// Indexed setters

void HDF5Matrix_integer_output_setRowIndexed_integer(void * ptr, size_t r, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::IntegerVector::iterator* out) {
    static_cast<HDF5IntOut*>(ptr)->set_row_indexed(r, n, *idx, *out);
    return;
}

void HDF5Matrix_integer_output_setColIndexed_integer(void * ptr, size_t c, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::IntegerVector::iterator* out) {
    static_cast<HDF5IntOut*>(ptr)->set_col_indexed(c, n, *idx, *out);
    return;
}

void HDF5Matrix_integer_output_setRowIndexed_numeric(void * ptr, size_t r, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::NumericVector::iterator* out) {
    static_cast<HDF5IntOut*>(ptr)->set_row_indexed(r, n, *idx, *out);
    return;
}

void HDF5Matrix_integer_output_setColIndexed_numeric(void * ptr, size_t c, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::NumericVector::iterator* out) {
    static_cast<HDF5IntOut*>(ptr)->set_col_indexed(c, n, *idx, *out);
    return;
}
