/**
 * @file IMatrixFile.hpp
 * @brief Interface for matrix files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_IMATRIXFILE_HPP
#define WILDRIVER_IMATRIXFILE_HPP




#include "IMatrixReader.hpp"
#include "IMatrixWriter.hpp"




namespace WildRiver
{


class IMatrixFile :
  public IMatrixReader,
  public IMatrixWriter
{
  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IMatrixFile() 
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
     * @brief Reset the current position in the matrix file to the first row.
     * Declared here to remove ambiguity.
     */
    virtual void firstRow() = 0;


    /**
     * @brief Get the filename/path of the current matrix.
     *
     * @return The filename/path. 
     */
    virtual std::string const & getFilename() const noexcept = 0;


};


}




#endif
