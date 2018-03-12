/**
 * @file VectorTextFile.cpp
 * @brief Implementation of abstract class for text based vector files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-07
 */




#include "VectorTextFile.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


VectorTextFile::VectorTextFile(
    std::string const & name) :
  TextFile(name)
{
  // do nothing
}


VectorTextFile::~VectorTextFile()
{
  // do nothing
}




}
