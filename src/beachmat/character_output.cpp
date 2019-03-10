#include "Rcpp.h"

#include "HDF5_writer.h"
#include "exports.h"
#include <vector>

template<>
char HDF5_writer<char, STRSXP>::get_empty() { return '\0'; }

template<>
Rcpp::RObject HDF5_writer<char, STRSXP>::get_firstval() {
    std::vector<char> first(default_type.getSize(), '\0');
    if (get_nrow() && get_ncol()) {
        extract_one(0, 0, first.data());
    }
    return Rcpp::StringVector::create(Rcpp::String(first.data()));
}

class HDF5_character_output {
public:
    HDF5_character_output(size_t nr, size_t nc, size_t strlen) : writer(nr, nc, strlen+1), bufsize(strlen+1),
        buffer(bufsize * std::max({ nr, nc, size_t(1) })) {}
    ~HDF5_character_output() = default;
    HDF5_character_output(const HDF5_character_output&) = default;
    HDF5_character_output& operator=(const HDF5_character_output&) = default;
    HDF5_character_output(HDF5_character_output&&) = default;
    HDF5_character_output& operator=(HDF5_character_output&&) = default;

    void get_row(size_t r, Rcpp::StringVector::iterator out, size_t first, size_t last) {
        char* ref=buffer.data();
        writer.extract_row(r, ref, first, last);
        for (size_t c=first; c<last; ++c, ref+=bufsize, ++out) {
            (*out)=ref;
        }
        return;
    }

    void get_col(size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) {
        char* ref=buffer.data();
        writer.extract_col(c, ref, first, last);
        for (size_t r=first; r<last; ++r, ref+=bufsize, ++out) {
            (*out)=ref;
        }
        return;
    }

    Rcpp::String get(size_t r, size_t c) {
        char* ref=buffer.data();
        writer.extract_one(r, c, ref);
        return ref;
    }

    void set_row(size_t r, Rcpp::StringVector::iterator in, size_t first, size_t last) {
        if (writer.get_ncol() + first >= last) { // ensure they can fit in 'buffer'; if not, it should trigger an error in insert_row().
            char* ref=buffer.data();
            for (size_t c=first; c<last; ++c, ref+=bufsize, ++in) {
                std::strncpy(ref, Rcpp::String(*in).get_cstring(), bufsize-1);
                ref[bufsize-1]='\0'; // strncpy only pads up to just before the last position.
            }
        }
        writer.insert_row(r, buffer.data(), first, last);
        return;
    }

    void set_col(size_t c, Rcpp::StringVector::iterator in, size_t first, size_t last) {
        if (writer.get_nrow() + first >= last) { // ensure they can fit in 'buffer'.
            char* ref=buffer.data();
            for (size_t r=first; r<last; ++r, ref+=bufsize, ++in) {
                std::strncpy(ref, Rcpp::String(*in).get_cstring(), bufsize-1);
                ref[bufsize-1]='\0';
            }
        }
        writer.insert_col(c, buffer.data(), first, last);
        return;
    }

    void set(size_t r, size_t c, Rcpp::String in) {
        char* ref=buffer.data();
        std::strncpy(ref, in.get_cstring(), bufsize-1);
        ref[bufsize-1]='\0';
        writer.insert_one(r, c, ref);
        return;
    }

    void set_col_indexed(size_t c, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::StringVector::iterator val) {
        if (buffer.size() < N*bufsize) {
            buffer.resize(N*bufsize);
        }

        char* ref=buffer.data();
        for (size_t i=0; i<N; ++i, ref+=bufsize, ++val) {
            std::strncpy(ref, Rcpp::String(*val).get_cstring(), bufsize-1);
            ref[bufsize-1]='\0';
        }

        writer.insert_col_indexed(c, N, idx, buffer.data());
        return;
    }

    void set_row_indexed(size_t r, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::StringVector::iterator val) {
        if (buffer.size() < N*bufsize) {
            buffer.resize(N*bufsize);
        }

        char* ref=buffer.data();
        for (size_t i=0; i<N; ++i, ref+=bufsize, ++val) {
            std::strncpy(ref, Rcpp::String(*val).get_cstring(), bufsize-1);
            ref[bufsize-1]='\0';
        }

        writer.insert_row_indexed(r, N, idx, buffer.data());
        return;
    }

    Rcpp::RObject yield() {
        return writer.yield();
    }
protected:
    HDF5_writer<char, STRSXP> writer;
    size_t bufsize;
    std::vector<char> buffer;
};

// Constructor, destructors and clones.

void * HDF5Matrix_character_output_create(size_t nr, size_t nc) {
    // TODO: allow users to adjust strlen by querying global variable.
    return static_cast<void*>(new HDF5_character_output(nr, nc, 10));
}

void HDF5Matrix_character_output_destroy(void * ptr) {
    delete static_cast<HDF5_character_output*>(ptr);
    return;
}

void * HDF5Matrix_character_output_clone(void * ptr) {
    HDF5_character_output* old=static_cast<HDF5_character_output*>(ptr);
    return static_cast<void*>(new HDF5_character_output(*old));
}

SEXP HDF5Matrix_character_output_yield(void * ptr) {
    return static_cast<HDF5_character_output*>(ptr)->yield();
}

// Basic getters

void HDF5Matrix_character_output_get(void * ptr, size_t r, size_t c, Rcpp::String* val) {
    *val=static_cast<HDF5_character_output*>(ptr)->get(r, c);
    return;
}

void HDF5Matrix_character_output_getRow(void * ptr, size_t r, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5_character_output*>(ptr)->get_row(r, *out, first, last);
    return;
}

void HDF5Matrix_character_output_getCol(void * ptr, size_t c, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5_character_output*>(ptr)->get_col(c, *out, first, last);
    return;
}

// Basic setters

void HDF5Matrix_character_output_set(void * ptr, size_t r, size_t c, Rcpp::String* val) {
    static_cast<HDF5_character_output*>(ptr)->set(r, c, *val);
    return;
}

void HDF5Matrix_character_output_setRow(void * ptr, size_t r, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5_character_output*>(ptr)->set_row(r, *out, first, last);
    return;
}

void HDF5Matrix_character_output_setCol(void * ptr, size_t c, Rcpp::StringVector::iterator* out, size_t first, size_t last) {
    static_cast<HDF5_character_output*>(ptr)->set_col(c, *out, first, last);
    return;
}

// Indexed setters

void HDF5Matrix_character_output_setRowIndexed(void * ptr, size_t r, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::StringVector::iterator* out) {
    static_cast<HDF5_character_output*>(ptr)->set_row_indexed(r, n, *idx, *out);
    return;
}

void HDF5Matrix_character_output_setColIndexed(void * ptr, size_t c, size_t n, Rcpp::IntegerVector::iterator* idx, Rcpp::StringVector::iterator* out) {
    static_cast<HDF5_character_output*>(ptr)->set_col_indexed(c, n, *idx, *out);
    return;
}


