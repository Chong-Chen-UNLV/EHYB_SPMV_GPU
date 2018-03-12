/**
 * @file VectorInHandle.hpp
 * @brief The VectorInHandle class.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2016
 * @version 1
 * @date 2016-10-25
 */




#ifndef WILDRIVER_VECTORINHANDLE_HPP
#define WILDRIVER_VECTORINHANDLE_HPP




#include <string>
#include <memory>
#include "IVectorReader.hpp"




namespace WildRiver 
{


class VectorInHandle
{
  public:
    /**
     * @brief Create a new file handle for reading vectors.
     *
     * @param name The filename/path of the file to read.
     */
    VectorInHandle(
        std::string const & name);


    /**
     * @brief Close the handle.
     */
    ~VectorInHandle();


    /**
     * @brief Get the size of the vector.
     *
     * @return The size of the vector.
     */
    ind_t getSize();


    /**
     * @brief Read the vector.
     *
     * @param vals The location to save the vector.
     * @param progress The variable to update as the vector is loaded from 0.0
     * to 1.0 (can be null).
     */
    void read(
        val_t * vals,
        double * progress = nullptr);


  private:
    std::unique_ptr<IVectorReader> m_reader;


    // disable copying
    VectorInHandle(
        VectorInHandle const & handle);
    VectorInHandle & operator=(
        VectorInHandle const & handle);




};




}




#endif
