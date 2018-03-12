/**
 * @file IGraphWriter.hpp
 * @brief Interface for writing graphs.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */





#ifndef WILDRIVER_IGRAPHWRITER_HPP
#define WILDRIVER_IGRAPHWRITER_HPP




#include <string>
#include <vector>


#include "base.h"
#include "MatrixEntry.hpp"
#include "Exception.hpp"




namespace WildRiver
{


class IGraphWriter
{
  private:


  protected:


  public:
    /**
     * @brief Destructor.
     */
    virtual ~IGraphWriter()
    {
    }


    /**
     * @brief Write a graph file from the given CSR structure.
     *
     * @param xadj The adjacency list pointer (length nvtxs+1).
     * @param adjncy The adjacency list (length nedges).
     * @param vwgt The vertex weights. 
     * @param adjwgt The edge weights.
     */
    virtual void write(
        ind_t const * xadj,
        dim_t const * adjncy,
        val_t const * vwgt,
        val_t const * adjwgt) = 0;


    /**
     * @brief Set the information for the graph. This must be called before
     * writing the graph.
     *
     * @param nvtxs The number of vertices.
     * @param nedges The number of edges (an undirected edge counts as two).
     * @param nvwgt The number of vertex weights (constraints).
     * @param ewgts Whether or not edge weights are present.
     */
    virtual void setInfo(
        dim_t nvtxs,
        ind_t nedges,
        int nvwgt,
        bool ewgts) = 0;


    /**
     * @brief Set the adjacency list and vertex weight of the next vertex.
     *
     * @param vwgts The vertex weights for this vertex.
     * @param list The adjacecny list.
     */
    virtual void setNextVertex(
        std::vector<val_t> const & vwgts,
        std::vector<matrix_entry_struct> const & list) = 0;


    /**
     * @brief Get the name of this filetype.
     *
     * @return The name of this filetype.
     */
    virtual std::string const & getName() const = 0;


};


}




#endif
