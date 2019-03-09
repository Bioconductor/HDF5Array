#ifndef BEACHMAT_HDF5_READER_H
#define BEACHMAT_HDF5_READER_H

#include "Rcpp.h"
#include "H5Cpp.h"

#include "beachmat/utils/dim_checker.h"
#include "beachmat/utils/utils.h"
#include "HDF5_utils.h"

#include <string>
#include <stdexcept>
#include <sstream>

/*** Class definition ***/

template<typename T, int RTYPE>
class HDF5_reader : public beachmat::dim_checker {
public:
    HDF5_reader(const Rcpp::RObject& incoming) : original(incoming),
        rowlist(H5::FileAccPropList::DEFAULT.getId()), collist(H5::FileAccPropList::DEFAULT.getId()) {

        auto classinfo=beachmat::get_class_package(incoming);
        if (classinfo.first!="HDF5Matrix" || classinfo.second!="HDF5Array") {
            throw std::runtime_error("matrix should be a HDF5Matrix object");
        }
    
        const Rcpp::RObject& h5_seed=beachmat::get_safe_slot(incoming, "seed");
        auto seedinfo=beachmat::get_class_package(h5_seed);
        if (seedinfo.first!="HDF5ArraySeed" || seedinfo.second!="HDF5Array") {
            throw std::runtime_error("'seed' slot in a HDF5Matrix object should be a HDF5ArraySeed object");
        }
    
        // Checking first value.
        const Rcpp::RObject& firstval=beachmat::get_safe_slot(h5_seed, "first_val");
        if (firstval.sexp_type()!=RTYPE) { 
            std::stringstream err;
            err << "'first_val' slot in a HDF5ArraySeed object should be " << beachmat::translate_type(RTYPE);
            throw std::runtime_error(err.str());
        }
    
        // Checking dimensions.
        if (!h5_seed.hasAttribute("dim")) { 
            throw std::runtime_error("HDF5ArraySeed object should have 'dim' attribute"); 
        }
        this->fill_dims(h5_seed.attr("dim"));
        const size_t& NC=this->ncol;
        const size_t& NR=this->nrow;
    
        // Checking names.
        try {
            filename=beachmat::make_to_string(beachmat::get_safe_slot(h5_seed, "filepath"));
        } catch (...) { 
            throw std::runtime_error("'filepath' slot in a HDF5ArraySeed object should be a string");
        }
        try {
            dataname=beachmat::make_to_string(beachmat::get_safe_slot(h5_seed, "name"));
        } catch (...) { 
            throw std::runtime_error("'name' slot in a HDF5ArraySeed object should be a string");
        }
        
        // Setting up the HDF5 accessors.
        hfile.openFile(filename.c_str(), H5F_ACC_RDONLY);
        hdata = hfile.openDataSet(dataname.c_str());
    
        H5::DataSpace hspace = hdata.getSpace();
        if (hspace.getSimpleExtentNdims()!=2) {
            throw std::runtime_error("data in HDF5 file is not a two-dimensional array");
        }
    
        hsize_t dims_out[2];
        hspace.getSimpleExtentDims(dims_out, NULL);
        if (dims_out[1]!=NR || dims_out[0]!=NC) { 
            throw std::runtime_error("dimensions in HDF5 file do not equal dimensions in the HDF5Matrix object");
        }
    
        hselect.set_dims(NR, NC);
    
        // Checking the type.
        auto curtype=hdata.getTypeClass();
        switch (RTYPE) {
            case REALSXP:
                if (curtype!=H5T_FLOAT) { 
                    throw std::runtime_error("data type in HDF5 file is not double");
                }
                break;
            case INTSXP: 
                if (curtype!=H5T_INTEGER) { 
                    throw std::runtime_error("data type in HDF5 file is not integer");
                }
                break;
            case LGLSXP:
                if (curtype!=H5T_INTEGER) { 
                    throw std::runtime_error("data type in HDF5 file is not logical");
                }
                break;
            case STRSXP:
                if (curtype!=H5T_STRING) { 
                    throw std::runtime_error("data type in HDF5 file is not character");
                }
                break;
        }

        // Setting the chunk cache parameters.
        calc_HDF5_chunk_cache_settings(this->nrow, this->ncol, hdata.getCreatePlist(), this->get_datatype(),
                onrow, oncol, rowokay, colokay, largerrow, largercol, rowlist, collist);
        return;
    }

    ~HDF5_reader() = default;
    HDF5_reader(const HDF5_reader&) = default;
    HDF5_reader& operator=(const HDF5_reader&) = default;
    HDF5_reader(HDF5_reader&&) = default;
    HDF5_reader& operator=(HDF5_reader&&) = default;

    template<typename X>
    void extract_row(size_t r, X* out, const H5::DataType& HDT, size_t first, size_t last) { 
        check_rowargs(r, first, last);
        reopen_HDF5_file_by_dim(filename, dataname, 
                hfile, hdata, H5F_ACC_RDONLY, rowlist, 
                onrow, oncol, largercol, rowokay);
        hselect.select_row(r, first, last);
        hdata.read(out, HDT, hselect.get_row_space(), hselect.get_mat_space());
        return;
    }

    template<typename X>
    void extract_col(size_t c, X* out, const H5::DataType& HDT, size_t first, size_t last) { 
        check_colargs(c, first, last);
        reopen_HDF5_file_by_dim(filename, dataname, 
                hfile, hdata, H5F_ACC_RDONLY, collist, 
                oncol, onrow, largerrow, colokay);
        hselect.select_col(c, first, last);
        hdata.read(out, HDT, hselect.get_col_space(), hselect.get_mat_space());
        return;
    }
        
    template<typename X>
    void extract_one(size_t r, size_t c, X* out, const H5::DataType& HDT) { 
        check_oneargs(r, c);
        hselect.select_one(r, c);
        hdata.read(out, HDT, hselect.get_one_space(), hselect.get_mat_space());
        return;
    }

    template<typename X>
    void extract_rows(Rcpp::IntegerVector::iterator rIt, size_t n, X* out, const H5::DataType& HDT, size_t first, size_t last) { 
        check_rowargs(0, first, last);
        check_row_indices(rIt, n);

        reopen_HDF5_file_by_dim(filename, dataname, 
                hfile, hdata, H5F_ACC_RDONLY, rowlist, 
                onrow, oncol, largercol, rowokay);

        hselect.clear_mat_space();
        for (size_t i=0; i<n; ++i, ++rIt) {
            hselect.select_row(*rIt, first, last, H5S_SELECT_OR);
        }

        hsize_t custom_dim[2];
        custom_dim[0]=last-first;
        custom_dim[1]=n;
        H5::DataSpace custom_space(2, custom_dim);
        hdata.read(out, HDT, custom_space, hselect.get_mat_space()); 
        return;
    }

    template<typename X>
    void extract_cols(Rcpp::IntegerVector::iterator cIt, size_t n, X* out, const H5::DataType& HDT, size_t first, size_t last) {
        check_colargs(0, first, last);
        check_col_indices(cIt, n);

        reopen_HDF5_file_by_dim(filename, dataname, 
                hfile, hdata, H5F_ACC_RDONLY, collist, 
                oncol, onrow, largerrow, colokay);

        hselect.clear_mat_space();
        for (size_t i=0; i<n; ++i, ++cIt) {
            hselect.select_col(*cIt, first, last, H5S_SELECT_OR);
        }

        hsize_t custom_dim[2];
        custom_dim[0]=n;
        custom_dim[1]=last-first;
        H5::DataSpace custom_space(2, custom_dim);
        hdata.read(out, HDT, custom_space, hselect.get_mat_space()); 
        return;
    }

    H5::DataType get_datatype() const {
        return hdata.getDataType();
    }
protected:
    Rcpp::RObject original;
    std::string filename, dataname;

    H5::H5File hfile;
    H5::DataSet hdata;
    HDF5_selector hselect;

    bool onrow, oncol;
    bool rowokay, colokay;
    bool largerrow, largercol;
    H5::FileAccPropList rowlist, collist;
};

#endif
