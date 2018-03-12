/**
 * @file IVectorWriter.hpp
 * @brief Interface for writing out vectors.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-07
 */




#ifndef WILDRIVER_IVECTORWRITER_HPP
#define WILDRIVER_IVECTORWRITER_HPP




#include <vector>

#include "base.h"




namespace WildRiver
{


class IVectorWriter
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IVectorWriter()
    {
    }


    /**
     * @brief Set the size of the vector.  
     *
     * @param size The new size of the vector.
     */
    virtual void setSize(
        ind_t size) = 0;


    /**
     * @brief Write the vector to the underlying medium. 
     *
     * @param vals The dense array of values in the vector.
     * @param progress The variable to update as teh vector is saved (can be
     * null).
     */
    virtual void write(
        val_t const * vals,
        double * progress) = 0;



};




}




#endif
