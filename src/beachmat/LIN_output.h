#ifndef BEACHMAT_LIN_OUTPUT_H
#define BEACHMAT_LIN_OUTPUT_H

#include "Rcpp.h"

#include "HDF5_writer.h"
#include "beachmat/utils/utils.h"

template<typename T, class V, int RTYPE>
class HDF5_lin_output {
private:
    HDF5_writer<T, RTYPE> writer;
public:
    HDF5_lin_output(size_t nr, size_t nc) : writer(nr, nc) {}
    ~HDF5_lin_output() = default;
    HDF5_lin_output(const HDF5_lin_output&) = default;
    HDF5_lin_output& operator=(const HDF5_lin_output&) = default;
    HDF5_lin_output(HDF5_lin_output&&) = default;
    HDF5_lin_output& operator=(HDF5_lin_output&&) = default;

    size_t get_nrow() const {
        return writer.get_nrow();
    }

    size_t get_ncol() const {
        return writer.get_ncol();
    }

    // Getters:
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

    T get(size_t r, size_t c) {
        T out;
        writer.extract_one(r, c, &out);
        return out;
    }

    void get_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
        writer.extract_row(r, &(*out), H5::PredType::NATIVE_INT32, first, last);
        return;
    }

    void get_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
        writer.extract_row(r, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
        return;
    }

    void get_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
        writer.extract_col(c, &(*out), H5::PredType::NATIVE_INT32, first, last);
        return;
    }

    void get_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
        writer.extract_col(c, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
        return;
    }

    // Setters:
    void set_col(size_t c, Rcpp::IntegerVector::iterator out) {
        set_col(c, out, 0, get_nrow());
        return;
    }

    void set_col(size_t c, Rcpp::NumericVector::iterator out) {
        set_col(c, out, 0, get_nrow());
        return;
    }

    void set_row(size_t r, Rcpp::IntegerVector::iterator out) {
        set_row(r, out, 0, get_ncol());
        return;
    }

    void set_row(size_t r, Rcpp::NumericVector::iterator out) {
        set_row(r, out, 0, get_ncol());
        return;
    }

    void set(size_t r, size_t c, T in) {
        writer.insert_one(r, c, &in);
        return;
    }

    void set_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
        writer.insert_row(r, &(*out), H5::PredType::NATIVE_INT32, first, last);
        return;
    }

    void set_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
        writer.insert_row(r, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
        return;
    }

    void set_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
        writer.insert_col(c, &(*out), H5::PredType::NATIVE_INT32, first, last);
        return;
    }

    void set_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
        writer.insert_col(c, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
        return;
    }

    void set_col_indexed(size_t c, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::IntegerVector::iterator val) {
        writer.insert_col_indexed(c, N, idx, val, H5::PredType::NATIVE_INT32);
        return;
    }
     
    void set_col_indexed(size_t c, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::NumericVector::iterator val) {
        writer.insert_col_indexed(c, N, idx, val, H5::PredType::NATIVE_DOUBLE);
        return;
    }

    void set_row_indexed(size_t r, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::IntegerVector::iterator val) {
        writer.insert_row_indexed(r, N, idx, val, H5::PredType::NATIVE_INT32);
        return;
    }

    void set_row_indexed(size_t r, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::NumericVector::iterator val) {
        writer.insert_row_indexed(r, N, idx, val, H5::PredType::NATIVE_DOUBLE);
        return;
    }

    // Other:
    Rcpp::RObject yield() {
        return writer.yield();
    }
};

#endif
