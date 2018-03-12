/**
 * @file VectorInHandle.cpp
 * @brief Implementation of the VectorInHandle class.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2016
 * @version 1
 * @date 2016-10-25
 */




#include "VectorInHandle.hpp"
#include "VectorFactory.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


VectorInHandle::VectorInHandle(
    std::string const & name) :
  m_reader(VectorFactory::make(name))
{
  // do nothing
}


VectorInHandle::~VectorInHandle()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


ind_t VectorInHandle::getSize()
{
  return m_reader->getSize();
}


void VectorInHandle::read(
    val_t * const vals,
    double * progress)
{
  m_reader->read(vals,progress);
}




}
