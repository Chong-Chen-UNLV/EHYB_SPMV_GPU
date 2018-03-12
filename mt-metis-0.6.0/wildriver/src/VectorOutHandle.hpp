/**
 * @file VectorOutHandle.hpp
 * @brief The VectorOutHandle class.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2016
 * @version 1
 * @date 2016-10-25
 */




#ifndef WILDRIVER_VECTOROUTHANDLE_HPP
#define WILDRIVER_VECTOROUTHANDLE_HPP




#include <string>
#include <memory>
#include "IVectorWriter.hpp"




namespace WildRiver
{


class VectorOutHandle
{
  public:
    /**
     * @brief Create a new file handle for writing vectors.
     *
     * @param name The filename/path of the file to read.
     */
    VectorOutHandle(
        std::string const & name);


    /**
     * @brief Close the handle;
     */
    ~VectorOutHandle();


    /**
     * @brief Write the vector.
     *
     * @param vals The location of the vector.
     * @param size The size of the vector.
     * @param progress A variable updated as the vector is written (from 0 to
     * 1.0). Can be null if no updates are required.
     */
    void write(
        val_t const * vals,
        ind_t size,
        double * progress = nullptr);


  private:
    std::unique_ptr<IVectorWriter> m_writer;


    // disable copying
    VectorOutHandle(
        VectorOutHandle const & handle);
    VectorOutHandle & operator=(
        VectorOutHandle const & handle);




};




}




#endif
