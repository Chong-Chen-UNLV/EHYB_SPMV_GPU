/**
 * @file GraphHandle.cpp
 * @brief Class for reading all matrix types. Uses PIMPL.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#include "GraphInHandle.hpp"
#include "GraphFactory.hpp"




namespace WildRiver
{


/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


GraphInHandle::GraphInHandle(
    std::string const & name) :
  m_reader(GraphFactory::make(name))
{
  // do nothing
}


GraphInHandle::~GraphInHandle()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void GraphInHandle::getInfo(
    dim_t & nvtxs,
    ind_t & nedges,
    int & nvwgt,
    bool & ewgts)
{
  m_reader->getInfo(nvtxs,nedges,nvwgt,ewgts);
}


void GraphInHandle::readGraph(
    ind_t * const xadj,
    dim_t * const adjncy,
    val_t * const vwgt,
    val_t * const adjwgt,
    double * progress)
{
  m_reader->read(xadj,adjncy,vwgt,adjwgt,progress);
}




}
