/**
 * @file GraphOutHandle.cpp
 * @brief Classfor writing all graph types.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#include "GraphOutHandle.hpp"
#include "GraphFactory.hpp"




namespace WildRiver
{



/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


GraphOutHandle::GraphOutHandle(
    std::string const & name) :
  m_writer(GraphFactory::make(name))
{
  // do nothing
}


GraphOutHandle::~GraphOutHandle()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void GraphOutHandle::writeGraph(
    ind_t const * const xadj,
    dim_t const * const adjncy,
    val_t const * const vwgt,
    val_t const * const adjwgt)
{
  m_writer->write(xadj,adjncy,vwgt,adjwgt);
}


void GraphOutHandle::setInfo(
    dim_t const nvtxs,
    ind_t const nedges,
    int const nvwgt,
    bool const ewgts)
{
  m_writer->setInfo(nvtxs,nedges,nvwgt,ewgts);
}


void GraphOutHandle::setNextVertex(
    std::vector<val_t> const & vwgts,
    std::vector<matrix_entry_struct> const & list)
{
  m_writer->setNextVertex(vwgts,list);
}




}
