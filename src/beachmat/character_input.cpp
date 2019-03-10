#include "Rcpp.h"

#include "HDF5_reader.h"
#include "beachmat/utils/dim_checker.h"
#include "exports.h"
#include <vector>

class HDF5_character_matrix {
public:
    HDF5_character_matrix(const Rcpp::RObject& incoming) : reader(incoming), str_type(reader.get_datatype()) {
        if (str_type.isVariableStr()) {
            throw std::runtime_error("variable-length strings not supported for HDF5_character_matrix");
        }

        bufsize=str_type.getSize();
        buffer.resize(bufsize * std::max({ reader.get_ncol(), reader.get_nrow(), size_t(1) }));
        return;
    }
    ~HDF5_character_matrix() = default;
    HDF5_character_matrix(const HDF5_character_matrix&) = default;
    HDF5_character_matrix& operator=(const HDF5_character_matrix&) = default;
    HDF5_character_matrix(HDF5_character_matrix&&) = default;
    HDF5_character_matrix& operator=(HDF5_character_matrix&&) = default;

    size_t get_nrow() const { return reader.get_nrow(); }
    size_t get_ncol() const { return reader.get_ncol(); }

    void get_row(size_t r, Rcpp::StringVector::iterator out, size_t first, size_t last) {
        char* ref=buffer.data();
        reader.extract_row(r, ref, str_type, first, last);
        for (size_t c=first; c<last; ++c, ref+=bufsize, ++out) {
            (*out)=ref;
        }
        return;
    }

    void get_col(size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) {
        char* ref=buffer.data();
        reader.extract_col(c, ref, str_type, first, last);
        for (size_t r=first; r<last; ++r, ref+=bufsize, ++out) {
            (*out)=ref;
        }
        return;
    }

    Rcpp::String get(size_t r, size_t c) {
        char* ref=buffer.data();
        reader.extract_one(r, c, ref, str_type);
        return ref;
    }

    void get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out, size_t first, size_t last) {
        const size_t required=bufsize * n * (last - first);
        if (required > buffer.size()) {
            beachmat::dim_checker::check_subset(first, last, reader.get_ncol(), "column");
            buffer.resize(required);
        }

        char* ref=buffer.data();
        reader.extract_rows(it, n, ref, str_type, first, last);
        for (size_t i=0; i<required; i+=bufsize, ref+=bufsize, ++out) {
            (*out)=ref;
        }
        return;
    }

    void get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out, size_t first, size_t last) {
        const size_t required=bufsize * n * (last - first);
        if (required > buffer.size()) {
            beachmat::dim_checker::check_subset(first, last, reader.get_nrow(), "row");
            buffer.resize(required);
        }

        char* ref=buffer.data();
        reader.extract_cols(it, n, ref, str_type, first, last);
        for (size_t i=0; i<required; i+=bufsize, ref+=bufsize, ++out) {
            (*out)=ref;
        }
        return;
    }

protected:
    HDF5_reader<char, STRSXP> reader;
    H5::DataType str_type;
    size_t bufsize;
    std::vector<char> buffer;
};

// Constructor, destructors and clones.

void * HDF5Matrix_character_input_create(SEXP incoming) {
    return static_cast<void*>(new HDF5_character_matrix(incoming));
}

void HDF5Matrix_character_input_destroy(void * ptr) {
    delete static_cast<HDF5_character_matrix*>(ptr);
    return;
}

void * HDF5Matrix_character_input_clone(void * ptr) {
    HDF5_character_matrix* old=static_cast<HDF5_character_matrix*>(ptr);
    return static_cast<void*>(new HDF5_character_matrix(*old));
}

// Basic getters

void HDF5Matrix_character_input_dim(void* ptr, size_t* nr, size_t* nc){
    HDF5_character_matrix* thing=static_cast<HDF5_character_matrix*>(ptr);
    *nr=thing->get_nrow();
    *nc=thing->get_ncol();
    return;
}

void HDF5Matrix_character_input_get(void * ptr, size_t r, size_t c, Rcpp::String* val) {
    *val=static_cast<HDF5_character_matrix*>(ptr)->get(r, c);
    return;
}

void HDF5Matrix_character_input_getRow(void * ptr, size_t r, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5_character_matrix*>(ptr)->get_row(r, *out, first, last);
    return;
}

void HDF5Matrix_character_input_getCol(void * ptr, size_t c, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5_character_matrix*>(ptr)->get_col(c, *out, first, last);
    return;
}

// Multi getters

void HDF5Matrix_character_input_getRows(void * ptr, Rcpp::IntegerVector::iterator* r, size_t n, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5_character_matrix*>(ptr)->get_rows(*r, n, *out, first, last);
    return;
}

void HDF5Matrix_character_input_getCols(void * ptr, Rcpp::IntegerVector::iterator* c, size_t n, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5_character_matrix*>(ptr)->get_cols(*c, n, *out, first, last);
    return;
}
