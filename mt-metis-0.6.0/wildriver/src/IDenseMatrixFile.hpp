/**
 * @file IDenseMatrixFile.hpp
 * @brief Interface for dense matrix files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-06
 */




#ifndef WILDRIVER_IDENSEMATRIXFILE_HPP
#define WILDRIVER_IDENSEMATRIXFILE_HPP




#include "IDenseMatrixReader.hpp"
#include "IDenseMatrixWriter.hpp"




namespace WildRiver
{


class IDenseMatrixFile :
  public IDenseMatrixReader,
  public IDenseMatrixWriter

{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IDenseMatrixFile()
    {
      // do nothing
    }


    /**
     * @brief Get the name of this matrix file type. Declared here to remove
     * ambiguity.
     *
     * @return The matrix file type name.
     */

    virtual std::string const & getName() const noexcept override = 0;


    /**
     * @brief Get the filename/path of the current matrix.
     *
     * @return The filename/path. 
     */
    virtual std::string const & getFilename() const noexcept = 0;




};




}




#endif
