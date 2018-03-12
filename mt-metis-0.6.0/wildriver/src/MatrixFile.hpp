/**
 * @file MatrixFile.hpp
 * @brief Base abstract class for matrix files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-05
 */




#ifndef WILDRIVER_MATRIXFILE_HPP
#define WILDRIVER_MATRIXFILE_HPP




#include <memory>
#include "Matrix.hpp"
#include "IMatrixFile.hpp"




namespace WildRiver
{


class MatrixFile :
  public Matrix,
  public IMatrixFile
{
  public:
    /**
     * @brief Default constructor which initializes a matrix with invalid
     * nrows, ncols, and nnz.
     */
    MatrixFile();


    /**
     * @brief Destructor.
     */
    virtual ~MatrixFile();


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
        ind_t & nnz) override;


    /**
     * @brief Set the number of rows, columns, and non-zeros in the martix.
     *
     * @param nrows The number of rows.
     * @param ncols The number of columns.
     * @param nnz The number of non-zeros.
     */
    virtual void setInfo(
        dim_t nrows,
        dim_t ncols,
        ind_t nnz) override;


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
     * @param progress The variable to update as the matrix is loaded (may be
     * null).
     */
    virtual void read(
        ind_t * rowptr,
        dim_t * rowind,
        val_t * rowval,
        double * progress) override = 0;


    /**
     * @brief Write the given CSR structure to teh file. The information for
     * the matrix must already be set.
     *
     * @param rowptr The row pointer indicating the start of each row.
     * @param rowind The row column indexs (i.e., for each element in a row,
     * the column index corresponding to that element).
     * @param rowval The row values.
     */
    virtual void write(
        ind_t const * rowptr,
        dim_t const * rowind,
        val_t const * rowval) override = 0;


  protected:
    /**
     * @brief Get the number of non-zero elements in this matrix.
     *
     * @return The number of non-zero elements in this matrix.
     */
    inline ind_t getNNZ() const noexcept
    {
      return m_nnz;
    }


    /**
     * @brief Set the number of non-zeroes in this matrix.
     *
     * @param n The new number of non-zeros.
     */
    inline void setNNZ(
        ind_t const n) noexcept
    {
      m_nnz = n;
    }


    /**
     * @brief Read the header of this matrix file. Populates internal fields
     * with the header information.
     */
    virtual void readHeader() = 0;


    /**
     * @brief Write the header of this matrix file. The header consists of
     * internal fields set by "setInfo()".
     */
    virtual void writeHeader() = 0; 


  private:
    /**
     * @brief The number of non-zeros in the matrix.
     */
    ind_t m_nnz;


    /**
     * @brief Indicating whether or not the matrix information has been set.
     */
    bool m_infoSet;





};




}




#endif
