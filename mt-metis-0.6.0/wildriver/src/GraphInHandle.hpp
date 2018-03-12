/**
 * @file GraphInHandle.hpp
 * @brief Class for reading all graph types. Uses PIMPL.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_GRAPHINHANDLE_HPP
#define WILDRIVER_GRAPHINHANDLE_HPP




#include <vector>
#include <memory>

#include "IGraphReader.hpp"




namespace WildRiver 
{


class GraphInHandle
{
  public:
    /**
     * @brief Create a new file handle for reading matrices.
     *
     * @param fname The filename/path of the file to read.
     */
    GraphInHandle(
        std::string const & fname);


    /**
     * @brief Close the handle.
     */
    ~GraphInHandle();


    /**
     * @brief Read the CSR structure of the graph.
     *
     * @param xadj The adjacency list pointer (length nvtxs+1).
     * @param adjncy The adjacency list (length nedges).
     * @param vwgt The vertex weights (length nvtxs*nvwgt). This may be NULL in
     * order to ignore vertex weights. If it is specified and the file does not
     * contain vertex weights, it will be filled with ones.
     * @param adjwgt The edge weights (length nedges). This may be NULL in
     * order to ignore edge weights. If it is specified and the file does not
     * contain edge weights, it will be filled with ones.
     * @param progress The variable to update as the graph is loaded (may be
     * null).
     */
    void readGraph(
        ind_t * xadj,
        dim_t * adjncy,
        val_t * vwgt,
        val_t * adjwgt,
        double * progress = nullptr);

  
    /**
     * @brief Get information about the graph.
     *
     * @param nvtxs The number of vertices.
     * @param nedges The number of edges (directed).
     * @param nvwgt The number of vertex weights (constraints).
     * @param ewgts Whether or not edge weights are specified.
     */
    void getInfo(
        dim_t & nvtxs,
        ind_t & nedges,
        int & nvwgt,
        bool & ewgts);


  private:
    /**
     * @brief A pointer to the underlying graph reader.
     */
    std::unique_ptr<IGraphReader> m_reader;


    /**
     * @brief Private copy constructor declared to disable copying.
     *
     * @param handle The handle to copy.
     */
    GraphInHandle(
        GraphInHandle const & handle);


    /**
     * @brief Private assignment operator declared to disable copying.
     *
     * @param handle The handle to copy.
     *
     * @return The new handle.
     */
    GraphInHandle & operator=(
        GraphInHandle const & handle);




};




}




#endif
