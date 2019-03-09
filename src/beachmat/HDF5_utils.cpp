#include "HDF5_utils.h"

#include <sstream>
#include <stdexcept>
#include <cmath>

/*******************************************
 * HDF5 functions, to be used for both
 * readers and writers.
 *******************************************/

/* HDF5 utilities. */

size_t get_cache_size_hard_limit () {
    return 2000000000;
}

/* This function computes the chunk cache settings for a HDF5 file
 * of a given dimension. It takes a bunch of HDF5_matrix/output 
 * members and modifies them by reference.
 */

void calc_HDF5_chunk_cache_settings (const size_t total_nrows, const size_t total_ncols, 
        const H5::DSetCreatPropList& cparms, const H5::DataType& default_type,
        bool& onrow, bool& oncol, bool& rowokay, bool& colokay, bool& largerrow, bool& largercol,
        H5::FileAccPropList& rowlist, H5::FileAccPropList& collist) {

    if (cparms.getLayout()!=H5D_CHUNKED) {
        // If contiguous, setting the flags to avoid reopening the file.
        onrow=true;
        oncol=true;
        rowokay=false;
        colokay=false;
        largerrow=false;
        largercol=false;
        return;
    }
    
    /* Setting up the chunk cache specification. */
    hsize_t chunk_dims[2];
    cparms.getChunk(2, chunk_dims);
    const size_t chunk_nrows=chunk_dims[1];
    const size_t chunk_ncols=chunk_dims[0];
    const size_t num_chunks_per_row=std::ceil(double(total_ncols)/chunk_ncols); // per row needs to divide by column dimensions.
    const size_t num_chunks_per_col=std::ceil(double(total_nrows)/chunk_nrows); 

    /* Everything is transposed, so hash indices are filled column-major. 
     * Here, we computing the lowest multiple of # chunks-per-col that is greater than # chunks-per-row, plus 1.
     * This ensures that two chunks in the same row/column do not have the same hash index.
     */
//  const size_t nslots = std::ceil(double(num_chunks_per_row)/num_chunks_per_col) * num_chunks_per_col + 1; 
    const size_t nslots = num_chunks_per_row * num_chunks_per_col; // fudge for https://forum.hdfgroup.org/t/unintended-behaviour-for-hash-values-during-chunk-caching/4869/.

    /* Computing the size of the cache required to store all chunks in each row or column.
     * The approach used below avoids overflow from computing eachchunk*num_Xchunks.
     */
    const size_t eachchunk=default_type.getSize() * chunk_nrows * chunk_ncols;
    const size_t nchunks_in_cache=get_cache_size_hard_limit()/eachchunk;
    rowokay=nchunks_in_cache >= num_chunks_per_row; 
    colokay=nchunks_in_cache >= num_chunks_per_col; 

    const size_t eachrow=eachchunk * num_chunks_per_row; 
    const size_t eachcol=eachchunk * num_chunks_per_col;
    largercol=eachcol >= eachrow;
    largerrow=eachrow >= eachcol;

    /* The first argument is ignored, according to https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html.
     * Setting w0 to 0 to evict the last used chunk; no need to worry about full vs partial reads here.
     */
    rowlist.setCache(10000, nslots, eachrow, 0);
    collist.setCache(10000, nslots, eachcol, 0);

    // File is not opened on either row or column yet.
    onrow=false; 
    oncol=false;
    return;
}

/* This function reopens the file with the chunk cache optimized for the
 * specified dimension. It needs a lot of internal components from the 
 * HDF5_matrix and HDF5_output objects to do so, including flags for the
 * other dimension (e.g., 'dim=>col' and 'other=>row' for column access).
 */

void reopen_HDF5_file_by_dim(const std::string& filename, const std::string& dataname, 
        H5::H5File& hfile, H5::DataSet& hdata, const unsigned& openmode, const H5::FileAccPropList& dimlist,
        bool& ondim, bool& onother, const bool& largerother, const bool& dimokay) {
    if (ondim || (onother && largerother)) {
        ; // Don't do anything, it's okay.
    } else if (!dimokay) {
        std::stringstream err;
        err << "cache size limit (" << get_cache_size_hard_limit() << ") exceeded for dim access, repack the file";
        throw std::runtime_error(err.str());
    } else {
        hdata.close();
        hfile.close();
        hfile.openFile(filename.c_str(), openmode, dimlist);
        hdata = hfile.openDataSet(dataname.c_str());
        ondim=true;
        onother=false;
    }
}

/***************************************************************************
 * This struct's methods set the rowspace and dataspace elements according to
 * the requested data access profile. We have column, row and single access.
 ***************************************************************************/

void HDF5_selector::set_dims(size_t NR, size_t NC) {
    h5_start[0]=0;
    h5_start[1]=0;

    one_count[0]=1;
    one_count[1]=1;
    one_space.setExtentSimple(1, one_count);
    one_space.selectAll();

    /* Using an 2D output space for 'row_space' and 'col_space', to match a 2D 'mat_space'.
     * Avoid unnecessary data rearrangements when going from 2D input to 1D output space.
     */
    row_count[0]=NC;
    row_count[1]=1;
    row_space.setExtentSimple(2, row_count);
    row_space.selectAll();

    col_count[0]=1;
    col_count[1]=NR;
    col_space.setExtentSimple(2, col_count);
    col_space.selectAll();

    hsize_t dims[2];
    dims[0]=NC;
    dims[1]=NR;
    mat_space.setExtentSimple(2, dims);
    return;
}

static const hsize_t zero_vals[]={0, 0};
const hsize_t* HDF5_selector::zero=zero_vals;

void HDF5_selector::select_row(size_t r, size_t start, size_t end, const H5S_seloper_t& op) { 
    hsize_t new_len = end - start;
    if (new_len != row_count[0]) { 
        row_count[0] = new_len;
        row_space.selectHyperslab(op, row_count, zero);
    }

    h5_start[0] = start;
    h5_start[1] = r;
    mat_space.selectHyperslab(op, row_count, h5_start);
    return;
}

void HDF5_selector::select_col(size_t c, size_t start, size_t end, const H5S_seloper_t& op) {
    hsize_t new_len = end - start;
    if (new_len != col_count[1]) { 
        col_count[1] = new_len;
        col_space.selectHyperslab(op, col_count, zero);
    }

    h5_start[0] = c;
    h5_start[1] = start;
    mat_space.selectHyperslab(op, col_count, h5_start);
    return;
}

void HDF5_selector::select_one(size_t r, size_t c, const H5S_seloper_t& op) {
    h5_start[0]=c;
    h5_start[1]=r;  
    mat_space.selectHyperslab(op, one_count, h5_start);
    return;
}

void HDF5_selector::select_indices(size_t n, const hsize_t* indices, const H5S_seloper_t& op) {
    mat_space.selectElements(op, n, indices);
    return;
}

void HDF5_selector::clear_mat_space () {
    mat_space.selectNone();
    return;
}

const H5::DataSpace& HDF5_selector::get_mat_space() const {
    return mat_space;
}

const H5::DataSpace& HDF5_selector::get_col_space() const {
    return col_space;
}

const H5::DataSpace& HDF5_selector::get_row_space() const {
    return row_space;
}

const H5::DataSpace& HDF5_selector::get_one_space() const {
    return one_space;
}

/* These overloaded functions return an output DataType for a 
 * given RTYPE. The second one also checks that the dataset is
 * of an appropriate class. Technically, I could put it in the 
 * input header, but it lives here with its overloaded sibling.
 */

H5::DataType set_HDF5_data_type (int RTYPE, size_t strlen) { 
    switch (RTYPE) {
        case REALSXP:
            return H5::PredType::NATIVE_DOUBLE;
        case INTSXP: 
            return H5::PredType::NATIVE_INT32;
        case LGLSXP:
            return H5::PredType::NATIVE_INT32;
        case STRSXP:
            return H5::StrType(0, strlen);
    }
    std::stringstream err;
    err << "unsupported sexptype '" << RTYPE << "'";
    throw std::runtime_error(err.str());
}
