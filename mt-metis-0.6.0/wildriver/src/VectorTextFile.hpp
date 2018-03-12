/**
 * @file VectorTextFile.hpp
 * @brief Abstract class for text based vector files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-07
 */




#ifndef WILDRIVER_VECTORTEXTFILE_HPP
#define WILDRIVER_VECTORTEXTFILE_HPP




#include "TextFile.hpp"
#include "VectorFile.hpp"




namespace WildRiver
{


class VectorTextFile :
  public TextFile,
  public VectorFile
{
  public:
    /**
     * @brief Open a vector file.
     *
     * @param name
     */
    VectorTextFile(
        std::string const & name);


    /**
     * @brief Virtual destructor.
     */
    virtual ~VectorTextFile();


    /**
     * @brief Get the filename/path of this file.
     *
     * @return The filename/path.
     */
    virtual std::string const & getFilename() const noexcept
    {
      return TextFile::getFilename();
    }




};




}




#endif
