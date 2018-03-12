/**
 * @file MatrixInHandle.hpp
 * @brief Class for reading all matrix types. Uses PIMPL.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_MATRIXINHANDLE_HPP
#define WILDRIVER_MATRIXINHANDLE_HPP




#include <memory>

#include "IMatrixReader.hpp"




namespace WildRiver 
{


class MatrixInHandle
{
  public:
    /**
     * @brief Create a new file handle for reading matrices.
     *
     * @param name The filename/path of the file to read.
     */
    MatrixInHandle(
        std::string const & name);


    /**
     * @brief Close the handle.
     */
    ~MatrixInHandle();


    /**
     * @brief Get the number of rows, columns, and non-zeros in the matrix.
     *
     * @param nrows The number of rows.
     * @param ncols The number of columns.
     * @param nnz THe number of non-zeros.
     */
    void getInfo(
        dim_t & nrows,
        dim_t & ncols,
        ind_t & nnz);


    /**
     * @brief Get the sparse matrix in CSR form. The pointers must be
     * pre-allocated to the sizes required by the info of the matrix 
     *
     * |rowptr| = nrows + 1
     * |rowind| = nnz
     * |rowval| = nnz
     *
     * @param rowptr The row pointer indicating the start of each row.
     * @param rowind The row column indexs (i.e., for each element in a row,
     * the column index corresponding to that element).
     * @param rowval The row values.
     * @param progress The variable to update as the matrix is loaded from 0.0
     * to 1.0 (can be null).
     */
    void readSparse(
        ind_t * rowptr,
        dim_t * rowind,
        val_t * rowval,
        double * progress = nullptr);


  private:
    std::unique_ptr<IMatrixReader> m_reader;


    // disable copying
    MatrixInHandle(
        MatrixInHandle const & handle);
    MatrixInHandle & operator=(
        MatrixInHandle const & handle);




};




}




#endif
