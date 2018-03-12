/**
 * @file MatrixHandle.cpp
 * @brief Class for reading all matrix types. Uses PIMPL.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#include "MatrixInHandle.hpp"
#include "MatrixFactory.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


MatrixInHandle::MatrixInHandle(
    std::string const & name) :
  m_reader(MatrixFactory::make(name))
  
{
  // do nothing
}


MatrixInHandle::~MatrixInHandle()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void MatrixInHandle::getInfo(
    dim_t & nrows,
    dim_t & ncols,
    ind_t & nnz)
{
  m_reader->getInfo(nrows,ncols,nnz);
}


void MatrixInHandle::readSparse(
    ind_t * rowptr,
    dim_t * rowind,
    val_t * rowval,
    double * progress)
{
  m_reader->read(rowptr,rowind,rowval,progress);
}




}
