/**
 * @file DenseMatrixFile.hpp
 * @brief Base abstract class for dense matrix files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-06
 */




#ifndef WILDRIVER_DENSEMATRIXFILE_HPP
#define WILDRIVER_DENSEMATRIXFILE_HPP




#include "IDenseMatrixFile.hpp"
#include "Matrix.hpp"




namespace WildRiver
{


class DenseMatrixFile :
  public Matrix,
  public IDenseMatrixFile
{
  public:
    /**
     * @brief Default constructor.
     */
    DenseMatrixFile();


    /**
     * @brief Cleanup associated memory.
     */
    ~DenseMatrixFile();


    /**
     * @brief Read the matrix into a dense data structure. The pointer must be
     * pre-allocated ot have room for all of the values of the matrix.
     *
     * |vals| = nrows*ncols
     *
     * On output, the values in the array will be in row major order. To access
     * the value in the ith row and jth column:
     *
     * a = vals[(i*ncols) + j];
     *
     * @param vals The dense array of values in the matrix.
     */
    virtual void read(
        val_t * vals) override = 0;


    /**
     * @brief Write the dense matrix out. The values in vals should be in row
     * major order. In other words, the value at the ith column and jth row
     * should be accessible via: 
     *
     * a = vals[(i*ncols) + j];
     *
     * @param vals The dense array of values in the matrix (of length
     * nrows*ncols).
     */
    virtual void write(
        val_t const * vals) override = 0;

};


}


#endif
