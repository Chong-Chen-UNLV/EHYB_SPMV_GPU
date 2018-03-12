/**
 * @file GraphFile.cpp
 * @brief Implementation of the GraphFile base abstract class.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-05
 */



#include "GraphFile.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


GraphFile::GraphFile() :
  m_ewgts(false),
  m_nvwgts(0)
{
  // do nothing
}


GraphFile::~GraphFile()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void GraphFile::getInfo(
    dim_t & nvtxs,
    ind_t & nedges,
    int & nvwgt,
    bool & ewgts)
{
  dim_t nrows, ncols;
  ind_t nnz;
  MatrixFile::getInfo(nrows,ncols,nnz);

  if (nrows != ncols) {
    throw BadParameterException(std::string("Graph files must be square: ") + \
        std::to_string(nrows) + std::string("x") + std::to_string(ncols));
  }

  nvtxs = nrows;
  nedges = nnz;

  nvwgt = getNumVertexWeights();
  ewgts = hasEdgeWeights();
}


void GraphFile::setInfo(
    dim_t const nvtxs,
    ind_t const nedges,
    int const nvwgt,
    bool const ewgts)
{
  setNumVertexWeights(nvwgt);
  hasEdgeWeights(ewgts);

  MatrixFile::setInfo(nvtxs,nvtxs,nedges);
}


void GraphFile::read(
    ind_t * rowptr,
    dim_t * rowind,
    val_t * rowval,
    double * progress)
{
  read(rowptr,rowind,NULL,rowval,progress);
}


void GraphFile::write(
    ind_t const * const rowptr,
    dim_t const * const rowind,
    val_t const * const rowval)
{
  // set weights
  if (rowval) {
    if (!hasEdgeWeights()) {
      hasEdgeWeights(true);

      // handle the case where we need to modify the header
      resetStream();
      writeHeader();
    }
  } else {
    hasEdgeWeights(false);
  }

  // make original call
  return write(rowptr,rowind,NULL,rowval);
}






}
