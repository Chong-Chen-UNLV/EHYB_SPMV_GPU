/**
 * @file IDenseMatrixWriter.hpp
 * @brief Interface for writing out dense matrices.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2016
 * @version 1
 * @date 2016-02-06
 */




#ifndef WILDRIVER_IDENSEMATRIXWRITER_HPP
#define WILDRIVER_IDENSEMATRIXWRITER_HPP




#include <string>
#include <vector>

#include "base.h"




namespace WildRiver
{


class IDenseMatrixWriter
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IDenseMatrixWriter()
    {
      // do nothing
    }


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
        val_t const * vals) = 0;


    /**
     * @brief Write the given row as the next row in the matrix.
     *
     * @param next The row to set.
     */
    virtual void setNextRow(
        std::vector<val_t> const & next) = 0;


    /**
     * @brief Get the name of this matrix writer type.
     *
     * @return The matrix writer type name.
     */
    virtual std::string const & getName() const noexcept = 0;




};




}




#endif
