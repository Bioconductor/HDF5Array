#ifndef BEACHMAT_HDF5_WRITER_H
#define BEACHMAT_HDF5_WRITER_H

#include "Rcpp.h"
#include "H5Cpp.h"

#include "beachmat/utils/dim_checker.h"
#include "HDF5_utils.h"

#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>

/*** Class definition ***/

template<typename T, int RTYPE>
class HDF5_writer : public beachmat::dim_checker {
public:
    HDF5_writer(size_t nr, size_t nc) : beachmat::dim_checker(nr, nc), 
            rowlist(H5::FileAccPropList::DEFAULT.getId()), collist(H5::FileAccPropList::DEFAULT.getId()) {

        // Pulling out settings.
        const Rcpp::Environment env=Rcpp::Environment::namespace_env("HDF5Array");
        Rcpp::Function fun=env["setupHDF5Matrix"];
        Rcpp::List collected=fun(Rcpp::IntegerVector::create(this->nrow, this->ncol), Rcpp::StringVector(beachmat::translate_type(RTYPE)));

        if (collected.size()!=4) { 
            throw std::runtime_error("output of setupHDF5Matrix should be a list of four elements");
        }
        fname=beachmat::make_to_string(collected[0]);
        dname=beachmat::make_to_string(collected[1]);
        Rcpp::IntegerVector r_chunks=collected[2];
        if (r_chunks.size()!=2) { 
            throw std::runtime_error("chunk dimensions should be an integer vector of length 2");
        }
        size_t chunk_nr=r_chunks[0];
        size_t chunk_nc=r_chunks[1];
        Rcpp::IntegerVector r_compress=collected[3];
        if (r_compress.size()!=1) {
            throw std::runtime_error("compression should be an integer scalar");
        }
        int compress=r_compress[0];

        // Opening the file, setting the type and creating the data set.
        hfile.openFile(fname.c_str(), H5F_ACC_RDWR);
        default_type=set_HDF5_data_type(RTYPE, 10); // TODO: allow users to set max string length.
        H5::DSetCreatPropList plist;
        const T empty=get_empty();
        plist.setFillValue(default_type, &empty);
            
        // Setting the chunk dimensions if not contiguous and has non-zero dimensions.
        if (compress>0 && nr && nc) {
            std::vector<hsize_t> chunk_dims(2);
            chunk_dims[0]=chunk_nc; // flipping them, as rhdf5 internally transposes it.
            chunk_dims[1]=chunk_nr;
            plist.setLayout(H5D_CHUNKED);
            plist.setChunk(2, chunk_dims.data());
            plist.setDeflate(compress);
        } else {
            plist.setLayout(H5D_CONTIGUOUS);
        }

        std::vector<hsize_t> dims(2);
        dims[0]=this->ncol; // Setting the dimensions (0 is column, 1 is row; internally transposed).
        dims[1]=this->nrow; 

        hselect.set_dims(this->nrow, this->ncol);
        hdata=hfile.createDataSet(dname.c_str(), default_type, hselect.get_mat_space(), plist); 

        // Setting logical attributes.
        if (RTYPE==LGLSXP) {
            std::string logi_attr="logical";
            H5::StrType str_type(0, logi_attr.size()+1);
            const hsize_t unit=1;
            H5::DataSpace att_space(1, &unit);
            H5::Attribute att = hdata.createAttribute("storage.mode", str_type, att_space);
            att.write(str_type, logi_attr);
        }

        // Setting the chunk cache parameters.
        calc_HDF5_chunk_cache_settings(this->nrow, this->ncol, hdata.getCreatePlist(), default_type, 
                onrow, oncol, rowokay, colokay, largerrow, largercol, rowlist, collist);
        return;
    }

    ~HDF5_writer() = default;
    HDF5_writer(const HDF5_writer&) = default;
    HDF5_writer& operator=(const HDF5_writer&) = default;
    HDF5_writer(HDF5_writer&&) = default;
    HDF5_writer& operator=(HDF5_writer&&) = default;

    // Setters:    
    template<typename X>
    void insert_col(size_t c, const X* in, const H5::DataType& HDT, size_t first, size_t last) {
        select_col(c, first, last);
        hdata.write(in, HDT, hselect.get_col_space(), hselect.get_mat_space());
        return;
    }
    
    void insert_col(size_t c, const T* in, size_t first, size_t last) {
        insert_col(c, in, default_type, first, last);
        return;
    }
    
    template<typename X>
    void insert_row(size_t c, const X* in, const H5::DataType& HDT, size_t first, size_t last) {
        select_row(c, first, last);
        hdata.write(in, HDT, hselect.get_row_space(), hselect.get_mat_space());
        return;
    }
    
    void insert_row(size_t c, const T* in, size_t first, size_t last) {
        insert_row(c, in, default_type, first, last);
        return;
    }
    
    void insert_one(size_t r, size_t c, T* in) {
        select_one(r, c);
        hdata.write(in, default_type, hselect.get_one_space(), hselect.get_mat_space());
        return;
    }

	// Indexed setters:
	void insert_col_indexed(size_t c, size_t n, const int* idx, const T* val) {
		insert_col_indexed(c, n, idx, val, default_type);    
	}

	template<typename X>
	void insert_col_indexed(size_t c, size_t n, const int* idx, const X* val, const H5::DataType& HDT) {
		if (!n) { return; }
		reopen_HDF5_file_by_dim(fname, dname,
				hfile, hdata, H5F_ACC_RDWR, collist, 
				oncol, onrow, largerrow, colokay);

		if (index_coords.size()/2 < n) {
			size_t N=std::max(n, this->nrow);
			index_coords.resize(N*2);
			hsize_t tmp_N=N;
			index_space.setExtentSimple(1, &tmp_N);
		}

		// Setting up the coordinate space.
		check_colargs(c);
		auto wsIt=index_coords.begin();
		for (size_t i=0; i<n; ++i, ++idx) {
			(*wsIt)=c;
			++wsIt;
			(*wsIt)=*idx;
			++wsIt;
		}

		hselect.select_indices(n, index_coords.data());

		// Performing a write operation.
		hsize_t tmp_count=n, tmp_start=0;
		index_space.selectHyperslab(H5S_SELECT_SET, &tmp_count, &tmp_start);
		hdata.write(val, HDT, index_space, hselect.get_mat_space());
		return;
	}

    void insert_row_indexed(size_t r, size_t n, const int* idx, const T* val) {
        insert_row_indexed(r, n, idx, val, default_type);    
    }

    template<typename X>
    void insert_row_indexed(size_t r, size_t n, const int* idx, const X* val, const H5::DataType& HDT) {
        if (!n) { return; }
        reopen_HDF5_file_by_dim(fname, dname, 
                hfile, hdata, H5F_ACC_RDWR, rowlist, 
                onrow, oncol, largercol, rowokay);

        if (index_coords.size()/2 < n) {
            size_t N=std::max(n, this->nrow);
            index_coords.resize(N*2);
            hsize_t tmp_N=N;
            index_space.setExtentSimple(1, &tmp_N);
        }

        // Setting up the coordinate space.
        check_rowargs(r);
        auto wsIt=index_coords.begin();
        for (size_t i=0; i<n; ++i, ++idx) {
            (*wsIt)=*idx;
            ++wsIt;
            (*wsIt)=r;
            ++wsIt;
        }

        hselect.select_indices(n, index_coords.data());

        // Performing a write operation.
        hsize_t tmp_count=n, tmp_start=0;
        index_space.selectHyperslab(H5S_SELECT_SET, &tmp_count, &tmp_start);
        hdata.write(val, HDT, index_space, hselect.get_mat_space());
        return;
    }

    // Getters:
    template<typename X>
    void extract_row(size_t r, X* out, const H5::DataType& HDT, size_t first, size_t last) { 
        select_row(r, first, last);
        hdata.read(out, HDT, hselect.get_row_space(), hselect.get_mat_space());
        return;
    } 

    void extract_row(size_t r, T* out, size_t first, size_t last) { 
        extract_row(r, out, default_type, first, last);
        return;
    } 

    template<typename X>
    void extract_col(size_t c, X* out, const H5::DataType& HDT, size_t first, size_t last) { 
        select_col(c, first, last);
        hdata.read(out, HDT, hselect.get_col_space(), hselect.get_mat_space());
        return;
    }

    void extract_col(size_t c, T* out, size_t first, size_t last) { 
        extract_col(c, out, default_type, first, last);
        return;
    }

    void extract_one(size_t r, size_t c, T* out) { 
        select_one(r, c);
        hdata.read(out, default_type, hselect.get_one_space(), hselect.get_mat_space());
        return;
    }

    // Returning a HDF5Matrix.
    Rcpp::RObject yield() {
        std::string seedclass="HDF5ArraySeed";
        Rcpp::S4 h5seed(seedclass);

        // Assigning to slots.
        if (!h5seed.hasSlot("filepath")) {
            throw std::runtime_error("missing 'filepath' slot in HDF5ArraySeed object");
        }
        h5seed.slot("filepath") = fname;

        if (!h5seed.hasSlot("name")) {
            throw std::runtime_error("missing 'name' slot in HDF5ArraySeed object");
        }
        h5seed.slot("name") = dname;

        if (!h5seed.hasSlot("dim")) {
            throw std::runtime_error("missing 'dim' slot in HDF5ArraySeed object");
        }
        h5seed.slot("dim") = Rcpp::IntegerVector::create(this->nrow, this->ncol);

        if (!h5seed.hasSlot("first_val")) {
            throw std::runtime_error("missing 'first_val' slot in HDF5ArraySeed object");
        }
        if (this->nrow && this->ncol) { 
            h5seed.slot("first_val") = get_firstval();
        } else {
            h5seed.slot("first_val") = Rcpp::Vector<RTYPE>(0); // empty vector.
        }

        // Assigning the seed to the HDF5Matrix.
        std::string matclass="HDF5Matrix";
        Rcpp::S4 h5mat(matclass);
        if (!h5mat.hasSlot("seed")) {
            throw std::runtime_error("missing 'seed' slot in HDF5Matrix object");
        }
        h5mat.slot("seed") = h5seed;
        return Rcpp::RObject(h5mat);
    }

protected:
    std::string fname, dname;

    H5::H5File hfile;
    H5::DataSet hdata;
    HDF5_selector hselect;

    H5::DataType default_type;
    void select_col(size_t c, size_t first, size_t last) {
        check_colargs(c, first, last);
        reopen_HDF5_file_by_dim(fname, dname,
                hfile, hdata, H5F_ACC_RDWR, collist, 
                oncol, onrow, largerrow, colokay);
        hselect.select_col(c, first, last);
        return;
    }

    void select_row(size_t r, size_t first, size_t last) {
        check_rowargs(r, first, last);
        reopen_HDF5_file_by_dim(fname, dname, 
                hfile, hdata, H5F_ACC_RDWR, rowlist, 
                onrow, oncol, largercol, rowokay);
        hselect.select_row(r, first, last);
        return;
    }

    void select_one(size_t r, size_t c) {
        check_oneargs(r, c);
        hselect.select_one(r, c);
        return;
    }

    bool onrow, oncol;
    bool rowokay, colokay;
    bool largerrow, largercol;
    H5::FileAccPropList rowlist, collist;

    // These functions are defined for each realized matrix type separately.
    static T get_empty();

    // In particular, this needs to be defined separately as strings are handled
    // very differently from LINs when creating any Rcpp::RObject.
    Rcpp::RObject get_firstval();

    // Objects for indexed access.
    std::vector<hsize_t> index_coords;
    H5::DataSpace index_space;   
};

#endif
