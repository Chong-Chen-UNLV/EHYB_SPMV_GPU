/**
 * @file MatrixEntry.hpp
 * @brief Structure a matrix entry.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_MATRIXENTRY_HPP
#define WILDRIVER_MATRIXENTRY_HPP




#include "base.h"




namespace WildRiver
{


/**
 * @brief Structure for holding the column index and value of an entry in a
 * sparse matrix.
 */
struct matrix_entry_struct 
{
  /**
   * @brief Column index of the entry.
   */
  dim_t ind;


  /**
   * @brief Numerical value of the entry.
   */
  val_t val;


};


}




#endif
