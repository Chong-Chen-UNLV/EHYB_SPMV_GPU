/**
 * @file IMatrixReader.hpp
 * @brief Interface for reading sparse matrices.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_IMATRIXREADER_HPP
#define WILDRIVER_IMATRIXREADER_HPP




#include <string>
#include <vector>

#include "base.h"
#include "MatrixEntry.hpp"
#include "Exception.hpp"




namespace WildRiver 
{


class IMatrixReader
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IMatrixReader() 
    {
    }


    /**
     * @brief Get the matrix in CSR form. The pointers must be
     * pre-allocated to the sizes required by the info of the matrix.
     *
     * |rowptr| = nrows + 1
     * |rowind| = nnz
     * |rowval| = nnz
     *
     * @param rowptr The row pointer indicating the start of each row.
     * @param rowind The row column indexs (i.e., for each element in a row,
     * the column index corresponding to that element).
     * @param rowval The row values.
     * @param progress The variable to update as the matrix is loaded (may be
     * null).
     */
    virtual void read(
        ind_t * rowptr,
        dim_t * rowind,
        val_t * rowval,
        double * progress) = 0;


    /**
     * @brief Get the number of rows, columns, and non-zeros in the matrix.
     *
     * @param nrows The number of rows.
     * @param ncols The number of columns.
     * @param nnz THe number of non-zeros.
     */
    virtual void getInfo(
        dim_t & nrows,
        dim_t & ncols,
        ind_t & nnz) = 0;


    /**
     * @brief Reset the current position in the matrix file to the first row.
     */
    virtual void firstRow() = 0;


    /**
     * @brief Get the next row in the matrix (adjacecny list in the graph).
     *
     * @param numNonZeros The number of non-zeros in the row (output).
     * @param columns The column of each non-zero entry (must be of length at
     * least the number of non-zero entries).
     * @param values The value of each non-zero entry (must be null or of 
     * length at least the number of non-zero entries).
     *
     * @return True if another row was found in the file.
     */
    virtual bool getNextRow(
        dim_t * numNonZeros,
        dim_t * columns,
        val_t * values) = 0;


    /**
     * @brief Get the name of this matrix file type.
     *
     * @return The matrix file type name.
     */
    virtual std::string const & getName() const noexcept = 0;


};




}




#endif
