/**
 * @file MatrixOutHandle.cpp
 * @brief Class for writing all matrix types. Uses PIMPL.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#include "MatrixOutHandle.hpp"
#include "MatrixFactory.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


MatrixOutHandle::MatrixOutHandle(
    std::string const & name) :
  m_writer(MatrixFactory::make(name))
{
  // do nothing
}


MatrixOutHandle::~MatrixOutHandle()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void MatrixOutHandle::setInfo(
    dim_t const nrows,
    dim_t const ncols,
    ind_t const nnz)
{
  m_writer->setInfo(nrows,ncols,nnz);
}


void MatrixOutHandle::writeSparse(
    ind_t const * const rowptr,
    dim_t const * const rowind,
    val_t const * const rowval)
{
  m_writer->write(rowptr,rowind,rowval);
}


void MatrixOutHandle::setNextRow(
    std::vector<matrix_entry_struct> const & next)
{
  m_writer->setNextRow(next);
}


}
