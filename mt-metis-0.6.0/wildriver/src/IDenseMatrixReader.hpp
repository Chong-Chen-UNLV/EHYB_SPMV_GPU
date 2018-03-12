/**
 * @file IDenseMatrixReader.hpp
 * @brief Interface for reading dense matrices.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-06
 */




#ifndef WILDRIVER_IDENSEMATRIXREADER_HPP
#define WILDRIVER_IDENSEMATRIXREADER_HPP




#include <vector>
#include <string>

#include "base.h"




namespace WildRiver
{


class IDenseMatrixReader
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IDenseMatrixReader()
    {
      // do nothing
    }


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
        val_t * vals) = 0;


    /**
     * @brief Read in the next row in the matrix.
     *
     * @param next The row to retrieve.
     *
     * @return True if the row was successfully read, and false if EOF is
     * reached.
     */
    virtual bool getNextRow(
        std::vector<val_t> & next) = 0;


    /**
     * @brief Get the name of this matrix file type.
     *
     * @return The matrix file type name.
     */
    virtual std::string const & getName() const noexcept = 0;


};




}




#endif
