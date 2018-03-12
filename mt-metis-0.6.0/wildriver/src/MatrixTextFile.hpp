/**
 * @file MatrixTextFile.hpp
 * @brief Abstract class for matrix files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_MATRIXTEXTFILE_HPP
#define WILDRIVER_MATRIXTEXTFILE_HPP




#include <vector>
#include <string>
#include <fstream>


#include "TextFile.hpp"
#include "MatrixFile.hpp"




namespace WildRiver 
{


class MatrixTextFile : 
    public TextFile,
    public MatrixFile 
{
  public:
    /**
     * @brief Open a new matrix file.
     *
     * @param fname The filename/path.
     */
    MatrixTextFile(
        std::string const & name);


    /**
     * @brief Virtual destructor.
     */
    virtual ~MatrixTextFile();


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
        double * progress) override;


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
        val_t const * rowval) override;


    /**
     * @brief Get the filename/path of this file.
     *
     * @return The filename/path.
     */
    virtual std::string const & getFilename() const noexcept
    {
      return TextFile::getFilename();
    }


  protected:
    /**
     * @brief Increase the current row number, and return the result.
     *
     * @return The current row number. 
     */
    dim_t incRow();


    /**
     * @brief Get the current row number.
     *
     * @return The current row number.
     */
    inline dim_t getCurrentRow() const noexcept
    {
      return current_row;
    }


    /**
     * @brief Set the current row number.
     *
     * @param row The new row number.
     */
    inline void setCurrentRow(
        dim_t const row) noexcept
    {
      current_row = row;
    }


  private:
    /**
     * @brief The number of non-zeros read in this file.
     */
    ind_t read_nnz;


    /**
     * @brief The row number pointed to by the current place in the filestream.
     */
    dim_t current_row;




};




}




#endif
