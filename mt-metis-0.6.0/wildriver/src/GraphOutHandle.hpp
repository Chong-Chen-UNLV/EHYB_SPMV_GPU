/**
 * @file GraphOutHandle.hpp
 * @brief Class for writing all graph types.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_GRAPHOUTHANDLE_HPP
#define WILDRIVER_GRAPHOUTHANDLE_HPP




#include <vector>
#include <memory>

#include "IGraphWriter.hpp"




namespace WildRiver
{


class GraphOutHandle
{
  public:
    /**
     * @brief Create a new file handle for writing matrices.
     *
     * @param fname The filename/path of the file to read.
     */
    GraphOutHandle(
        std::string const & fname);


    /**
     * @brief Close the matrix handle.
     */
    ~GraphOutHandle();


    /**
     * @brief Write a graph file from the given CSR structure.
     *
     * @param xadj The adjacency list pointer (length nvtxs+1).
     * @param adjncy The adjacency list (length nedges).
     * @param vwgt The vertex weights. 
     * @param adjwgt The edge weights.
     */
    void writeGraph(
        ind_t const * xadj,
        dim_t const * adjncy,
        val_t const * vwgt,
        val_t const * adjwgt);


    /**
     * @brief Set the information for the graph. This must be called before
     * writing the graph.
     *
     * @param nvtxs The number of vertices.
     * @param nedges The number of edges (an undirected edge counts as two).
     * @param nvwgt The number of vertex weights (constraints).
     * @param ewgts Whether or not edge weights are present.
     */
    void setInfo(
        dim_t nvtxs,
        ind_t nedges,
        int nvwgt,
        bool ewgts);


    /**
     * @brief Set the adjacency list and vertex weight of the next vertex.
     *
     * @param vwgts The vertex weights for this vertex.
     * @param list The adjacecny list.
     */
    void setNextVertex(
        std::vector<val_t> const & vwgts,
        std::vector<matrix_entry_struct> const & list);


  private:
    /**
     * @brief A pointer to the underlying graph writer.
     */
    std::unique_ptr<IGraphWriter> m_writer;


    // disable copying
    GraphOutHandle(
        GraphOutHandle const & handle);
    GraphOutHandle & operator=(
        GraphOutHandle const & handle);




};


}




#endif
