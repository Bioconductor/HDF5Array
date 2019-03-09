#ifndef BEACHMAT_LIN_MATRIX_H
#define BEACHMAT_LIN_MATRIX_H

#include "Rcpp.h"

#include "HDF5_reader.h"
#include "beachmat/utils/utils.h"

/* HDF5Matrix of LINs */

template<typename T, int RTYPE>
class HDF5_lin_matrix { 
private:
    HDF5_reader<T, RTYPE> reader;
public:
    HDF5_lin_matrix(const Rcpp::RObject& incoming) : reader(incoming) {}

    ~HDF5_lin_matrix() = default;
    HDF5_lin_matrix(const HDF5_lin_matrix&) = default;
    HDF5_lin_matrix& operator=(const HDF5_lin_matrix&) = default;
    HDF5_lin_matrix(HDF5_lin_matrix&&) = default;
    HDF5_lin_matrix& operator=(HDF5_lin_matrix&&) = default;

    size_t get_nrow() const {
        return reader.get_nrow();
    }

    size_t get_ncol() const {
        return reader.get_ncol();
    }

    T get(size_t r, size_t c) {
        T val;
        reader.extract_one(r, c, &val, reader.get_datatype());
        return val;
    }

    // Basic getters.
    void get_col(size_t c, Rcpp::IntegerVector::iterator out) {
        get_col(c, out, 0, get_nrow());
        return;
    }

    void get_col(size_t c, Rcpp::NumericVector::iterator out) {
        get_col(c, out, 0, get_nrow());
        return;
    }

    void get_row(size_t r, Rcpp::IntegerVector::iterator out) {
        get_row(r, out, 0, get_ncol());
        return;
    }

    void get_row(size_t r, Rcpp::NumericVector::iterator out) {
        get_row(r, out, 0, get_ncol());
        return;
    }

    void get_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
        reader.extract_col(c, static_cast<int*>(out), H5::PredType::NATIVE_INT32, first, last);
        return;
    }

    void get_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
        reader.extract_col(c, static_cast<double*>(out), H5::PredType::NATIVE_DOUBLE, first, last);
        return;
    }

    void get_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
        reader.extract_row(r, static_cast<int*>(out), H5::PredType::NATIVE_INT32, first, last);
        return;
    }

    void get_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
        reader.extract_row(r, static_cast<double*>(out), H5::PredType::NATIVE_DOUBLE, first, last);
        return;
    }

    // Multi getters.
    void get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out) {
        get_cols(it, n, out, 0, get_nrow());
        return;
    }

    void get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out) {
        get_cols(it, n, out, 0, get_nrow());
        return;
    }

    void get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out) {
        get_rows(it, n, out, 0, get_ncol());
        return;
    }

    void get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out) {
        get_rows(it, n, out, 0, get_ncol());
        return;
    }

    void get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
        reader.extract_cols(it, n, static_cast<int*>(out), H5::PredType::NATIVE_INT32, first, last);
        return;
    }

    void get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
        reader.extract_cols(it, n, static_cast<double*>(out), H5::PredType::NATIVE_DOUBLE, first, last);
        return;
    }

    void get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
        reader.extract_rows(it, n, static_cast<int*>(out), H5::PredType::NATIVE_INT32, first, last);
        return;
    }

    void get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
        reader.extract_rows(it, n, static_cast<double*>(out), H5::PredType::NATIVE_DOUBLE, first, last);
        return;
    }
};

#endif
