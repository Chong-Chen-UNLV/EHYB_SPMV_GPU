/**
 * @file MatrixFile.cpp
 * @brief Implementation of the base abstract class for matrix files. 
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-05
 */




#include "MatrixFile.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


MatrixFile::MatrixFile() :
  m_nnz(NULL_IND),
  m_infoSet(false)
{
  // do nothing
}


MatrixFile::~MatrixFile()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void MatrixFile::getInfo(
    dim_t & nrows,
    dim_t & ncols,
    ind_t & nnz)
{
  // see if need to read the header
  if (!m_infoSet) {
    readHeader();
  }

  // set values
  nrows = getNumRows();
  ncols = getNumCols();
  nnz = getNNZ();

  m_infoSet = true;
}


void MatrixFile::setInfo(
    dim_t const nrows,
    dim_t const ncols,
    ind_t const nnz)
{
  setNumRows(nrows);
  setNumCols(ncols);
  setNNZ(nnz);

  m_infoSet = true;

  writeHeader();
}




}
