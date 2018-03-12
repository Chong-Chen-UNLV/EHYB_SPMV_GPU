/**
 * @file GraphTextFile.cpp
 * @brief Abstract class for reading and writing graphs.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 */




#include "GraphTextFile.hpp"
#include "MetisFile.hpp"
#include "CSRFile.hpp"
#include "MatrixGraphFile.hpp"




namespace WildRiver
{


/******************************************************************************
* PROTECTED FUNCTIONS *********************************************************
******************************************************************************/


dim_t GraphTextFile::incVertex()
{
  if (m_currentVertex >= getNumVertices()) {
    throw BadFileStateException(std::string("Attempt to increase current " \
        "vertex beyond the number of vertices: ") + \
        std::to_string(getNumRows()));
  }

  ++m_currentVertex;

  return m_currentVertex;
}




/******************************************************************************
* CONSTRUCTOR / DESTRUCTOR ****************************************************
******************************************************************************/


GraphTextFile::GraphTextFile(
    std::string const & name) :
  TextFile(name),
  m_currentVertex(0)
{
  // do nothing
}


GraphTextFile::~GraphTextFile()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void GraphTextFile::read(
    ind_t * const xadj,
    dim_t * const adjncy,
    val_t * const vwgt,
    val_t * const adjwgt,
    double * progress)
{
  dim_t const nvtxs = getNumVertices();
  dim_t const ncon = getNumVertexWeights();

  dim_t const interval = nvtxs > 100 ? nvtxs / 100 : 1;
  double const increment = 1/100.0;

  xadj[0] = 0;
  for (dim_t i=0; i<nvtxs; ++i) {
    dim_t degree;

    val_t * const vwgtStart = (ncon > 0 && vwgt) ? vwgt+(i*ncon) : nullptr;
    dim_t * const adjncyStart = adjncy+xadj[i];
    val_t * const adjwgtStart = adjwgt ? adjwgt+xadj[i] : nullptr;

    // retrieve row
    if (!getNextVertex(vwgtStart,&degree,adjncyStart,adjwgtStart)) {
      throw BadFileException(std::string("Premature end of file: ") + \
          std::to_string(i) + std::string("/") + std::to_string(nvtxs) + \
          std::string(" vertices found."));
    }

    // handle vertex weights
    if (ncon == 0 && vwgt) {
      // set unit vertex weights
      vwgt[i] = 1;
    }

    xadj[i+1] = xadj[i]+degree;

    if (progress != nullptr && i % interval == 0) {
      *progress += increment;
    }
  }
}


void GraphTextFile::write(
    ind_t const * const xadj,
    dim_t const * const adjncy,
    val_t const * const vwgt,
    val_t const * const adjwgt)
{
  dim_t const nvtxs = getNumVertices();
  dim_t const ncon = getNumVertexWeights();
  bool const ewgts = hasEdgeWeights();

  for (dim_t i=0; i<nvtxs; ++i) {
    std::vector<val_t> vwgts;
    std::vector<matrix_entry_struct> list;

    // build vertex weight vector
    for (dim_t k=0; k<ncon; ++k) {
      vwgts.push_back(vwgt[(i*ncon)+k]);
    }

    // build edges
    for (ind_t j=xadj[i]; j<xadj[i+1]; ++j) {
      matrix_entry_struct e;
      e.ind = adjncy[j];
      if (ewgts) {
        if (adjwgt) {
          e.val = adjwgt[j];
        } else {
          e.val = 1;
        }
      }
      list.push_back(e);
    }

    // set the vertex
    setNextVertex(vwgts,list);
  }
}




}
