/**
 * @file GraphFactory.cpp
 * @brief Implementation of GraphFactory for instantiating graph files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-05
 */



#include "GraphFactory.hpp"
#include "MetisFile.hpp"
#include "CSRFile.hpp"
#include "MatrixGraphFile.hpp"
#include "MatrixFactory.hpp"




namespace WildRiver
{


/******************************************************************************
* STATIC FUNCTIONS ************************************************************
******************************************************************************/


std::unique_ptr<IGraphFile> GraphFactory::make(
    std::string const & name)
{
  std::unique_ptr<IGraphFile> file;
  // determine what type of reader to instantiate based on extension
  if (MetisFile::hasExtension(name)) {
    file.reset(new MetisFile(name));
  } else {
    // need to wrap it with an adapter
    std::unique_ptr<IMatrixFile> matPtr(MatrixFactory::make(name));
    file.reset(new MatrixGraphFile(matPtr));
  }
 
  return file;
}




}
