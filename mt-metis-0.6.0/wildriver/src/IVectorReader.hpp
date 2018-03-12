/**
 * @file IVectorReader.hpp
 * @brief Interface for reading in vectors.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-07
 */




#ifndef WILDRIVER_IVECTORREADER_HPP
#define WILDRIVER_IVECTORREADER_HPP




#include "base.h"




namespace WildRiver
{


class IVectorReader
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IVectorReader()
    {
      // do nothing
    }


    /**
     * @brief Get the size of the vector. May alter the internal state of the
     * reader.
     *
     * @return The size of the vector.
     */
    virtual ind_t getSize() = 0;


   /**
     * @brief Read the values of the vector. 
     *
     * @param vals The values in the vector (output).
     * @param progress The variable to update as teh vector is loaded (can be
     * null).
     */
    virtual void read(
        val_t * vals,
        double * progress) = 0;




};




}




#endif
