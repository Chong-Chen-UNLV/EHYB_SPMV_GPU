/**
 * @file MatrixFactory.cpp
 * @brief Implementation of the MatrixFactory class.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-05
 */




#include "MatrixFactory.hpp"
#include "CSRFile.hpp"
#include "MetisFile.hpp"




namespace WildRiver
{


/******************************************************************************
* STATIC FUNCTIONS ************************************************************
******************************************************************************/


std::unique_ptr<IMatrixFile> MatrixFactory::make(
    std::string const & name)
{
  std::unique_ptr<IMatrixFile> file;

  // determine what type of reader to instantiate based on extension
  if (MetisFile::hasExtension(name)) {
    file.reset(new MetisFile(name));
  } else if (CSRFile::hasExtension(name)) {
    file.reset(new CSRFile(name));
  } else {
    throw UnknownExtensionException(std::string("Unknown filetype: ") + name);
  }
 
  return file;
}




}
