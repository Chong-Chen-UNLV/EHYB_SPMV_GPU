/**
 * @file IVectorFile.hpp
 * @brief Interface for vector files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-07
 */




#ifndef WILDIRVER_IVECTORFILE_HPP
#define WILDIRVER_IVECTORFILE_HPP




#include "IVectorReader.hpp"
#include "IVectorWriter.hpp"




namespace WildRiver
{


class IVectorFile :
  public IVectorReader,
  public IVectorWriter
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IVectorFile()
    {
    }


    /**
     * @brief Get the filename/path of the current vector file.
     *
     * @return The filename/path. 
     */
    virtual std::string const & getFilename() const noexcept = 0;


};




}




#endif
