/**
 * @file MatrixOutHandle.hpp
 * @brief Class for writing all matrix types. Uses PIMPL.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_MATRIXOUTHANDLE_HPP
#define WILDRIVER_MATRIXOUTHANDLE_HPP




#include <vector>
#include <memory>

#include "IMatrixWriter.hpp"




namespace WildRiver 
{


class MatrixOutHandle
{
  public:
    /**
     * @brief Create a new handle for writing out matrices.
     *
     * @param name The filename to write the matrix to.
     */
    MatrixOutHandle(
        std::string const & name);


    /**
     * @brief Destructor.
     */
    ~MatrixOutHandle();


    /**
     * @brief Set information about the matrix. This function must be called
     * before writeSparse().
     *
     * @param numRows The number of rows in the matrix.
     * @param numCols The number of columns in the matrix.
     * @param numNonZeros The number of non-zeros in the matrix. 
     */
    void setInfo(
        dim_t nrows,
        dim_t ncols,
        ind_t nnz);


    /**
     * @brief Write the data of teh matrix. Thi function must be called after
     * setInfo().
     *
     * @param rowptr The beginning index of each row.
     * @param rowind The column each entry is located in.
     * @param rowval The value of each entry.
     */
    void writeSparse(
        ind_t const * rowptr,
        dim_t const * rowind,
        val_t const * rowval);


    /**
     * @brief Set the entries in the next row of the matrix.
     *
     * @param next The entries.
     *
     * @return 
     */
    void setNextRow(
        std::vector<matrix_entry_struct> const & next);


  private:
    std::unique_ptr<IMatrixWriter> m_writer;

    // disable copying
    MatrixOutHandle(
        MatrixOutHandle const & handle);
    MatrixOutHandle & operator=(
        MatrixOutHandle const & handle);





};




}




#endif
