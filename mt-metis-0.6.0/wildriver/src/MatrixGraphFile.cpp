/**
 * @file MatrixGraphFile.cpp
 * @brief An adapter class for treating matrices as graphs.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#include "MatrixGraphFile.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


MatrixGraphFile::MatrixGraphFile(
    std::unique_ptr<IMatrixFile>& file) :
  m_file(std::move(file))
{
  // do nothing
}


MatrixGraphFile::~MatrixGraphFile()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/

void MatrixGraphFile::getInfo(
    dim_t & nvtxs,
    ind_t & nedges,
    int & nvwgt,
    bool & ewgts)
{
  dim_t nrows, ncols;
  ind_t nnz;

  m_file->getInfo(nrows,ncols,nnz);

  nvtxs = nrows;
  nedges = nnz;
  
  // matrices never have vertex weights and always have edge weights
  nvwgt = 0;
  ewgts = true;
}


void MatrixGraphFile::setInfo(
    dim_t const nvtxs,
    ind_t const nedges,
    int const nvwgt,
    bool const ewgts)
{
  if (nvwgt != 0) {
    throw BadParameterException("MatrixGraphFile cannot handle vertex " \
        "weights");
  }

  // TODO: add an internal state for edge weights, such that they can be
  // surpessed via this function.

  m_file->setInfo(nvtxs,nvtxs,nedges);
}


void MatrixGraphFile::read(
    ind_t * const xadj,
    dim_t * const adjncy,
    val_t * const vwgt,
    val_t * const adjwgt,
    double * progress)
{
  dim_t nrows, ncols;
  ind_t nnz;

  m_file->getInfo(nrows,ncols,nnz);

  // fill vertex weights with ones if requested
  if (vwgt) {
    for (dim_t i=0; i<nrows; ++i) {
      vwgt[i] = 1;
    }
  }

  m_file->read(xadj,adjncy,adjwgt,progress);
}


void MatrixGraphFile::write(
    ind_t const * const xadj,
    dim_t const * const adjncy,
    val_t const * const vwgt,
    val_t const * const adjwgt)
{
  m_file->write(xadj,adjncy,adjwgt);
}


void MatrixGraphFile::firstVertex()
{
  m_file->firstRow();
}


bool MatrixGraphFile::getNextVertex(
    val_t *,
    dim_t * const degree,
    dim_t * const edgeDests,
    val_t * const edgeWeights)
{
  return m_file->getNextRow(degree,edgeDests,edgeWeights);
}


void MatrixGraphFile::setNextVertex(
    std::vector<val_t> const & vwgts,
    std::vector<matrix_entry_struct> const & list)
{
  m_file->setNextRow(list);
}


std::string const & MatrixGraphFile::getName() const noexcept
{
  return m_file->getName();
}




}
