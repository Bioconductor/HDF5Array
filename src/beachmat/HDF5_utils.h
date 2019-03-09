#ifndef BEACHMAT_HDF5_UTILS_H
#define BEACHMAT_HDF5_UTILS_H

#include "H5Cpp.h"

#include "beachmat/utils/utils.h"

#include <string>

size_t get_cache_size_hard_limit();

void calc_HDF5_chunk_cache_settings (const size_t, const size_t, const H5::DSetCreatPropList&, const H5::DataType&,
        bool&, bool&, bool&, bool&, bool&, bool&,
        H5::FileAccPropList&, H5::FileAccPropList&);

void reopen_HDF5_file_by_dim(const std::string&, const std::string&, 
        H5::H5File&, H5::DataSet&, const unsigned&, const H5::FileAccPropList&,
        bool&, bool&, const bool&, const bool&);

struct HDF5_selector {
    void set_dims(size_t, size_t);   
    void select_row(size_t, size_t, size_t, const H5S_seloper_t& op=H5S_SELECT_SET);
    void select_col(size_t, size_t, size_t, const H5S_seloper_t& op=H5S_SELECT_SET);
    void select_one(size_t, size_t, const H5S_seloper_t& op=H5S_SELECT_SET);
    void select_indices(size_t, const hsize_t*, const H5S_seloper_t& op=H5S_SELECT_SET);

    void clear_mat_space();

    const H5::DataSpace& get_mat_space() const;
    const H5::DataSpace& get_col_space() const;
    const H5::DataSpace& get_row_space() const;
    const H5::DataSpace& get_one_space() const;
private:
    H5::DataSpace col_space, row_space, one_space, mat_space;
    hsize_t h5_start[2], col_count[2], row_count[2], one_count[2];
    static const hsize_t* zero;
};

H5::DataType set_HDF5_data_type (int, size_t);

#endif
