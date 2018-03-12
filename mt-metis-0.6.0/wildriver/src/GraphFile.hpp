/**
 * @file GraphFile.hpp
 * @brief Base abstract class for graph files. 
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-05
 */




#ifndef WILDRIVER_GRAPHFILE_HPP
#define WILDRIVER_GRAPHFILE_HPP




#include "IGraphFile.hpp"
#include "MatrixFile.hpp"




namespace WildRiver
{


class GraphFile :
  public MatrixFile,
  public IGraphFile
{
  public:
    /**
     * @brief Default constructor which initializes a graph with an invalid
     * number of vertices and edges.
     */
    GraphFile();


    /**
     * @brief Destructor.
     */
    virtual ~GraphFile();


    /**
     * @brief Get information about the graph.
     *
     * @param nvtxs The number of vertices.
     * @param nedges The number of edges (directed).
     * @param nvwgt The number of vertex weights (constraints).
     * @param ewgts Whether or not edge weights are specified.
     */
    virtual void getInfo(
        dim_t & nvtxs,
        ind_t & nedges,
        int & nvwgt,
        bool & ewgts) override;


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
        bool ewgts) override;


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
    virtual void read(
        ind_t * xadj,
        dim_t * adjncy,
        val_t * vwgt,
        val_t * adjwgt,
        double * progress) override = 0;


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
        val_t const * adjwgt) override = 0;


    /**
     * @brief Get the graph in CSR form. The pointers must be
     * pre-allocated to the sizes required by the info of the matrix 
     *
     * |rowptr| = nrows + 1
     * |rowind| = nnz
     * |rowval| = nnz
     *
     * @param rowptr The row pointer indicating the start of each row.
     * @param rowind The row column indexs (i.e., for each element in a row,
     * the column index corresponding to that element).
     * @param rowval The row values.
     * @param progress The variable to update as the matrix is loaded (may be
     * null).
     */
    virtual void read(
        ind_t * rowptr,
        dim_t * rowind,
        val_t * rowval,
        double * progress) override;


    /**
     * @brief Write the CSR structure to the graph file.
     *
     * @param rowptr The row pointer / xadj.
     * @param rowind The row index / adjncy.
     * @param rowval The row val / adjwgt.
     */
    virtual void write(
        ind_t const * rowptr,
        dim_t const * rowind,
        val_t const * rowval) override;


   protected:
    /**
     * @brief Get the number of vertices in this graph.
     *
     * @return The number of vertices in this graph.
     */
    inline dim_t getNumVertices() const noexcept
    {
      return getNumRows();
    }


    /**
     * @brief Set the number of vertices in this graph.
     *
     * @param nvtxs The number of vertices in this graph.
     */
    inline void setNumVertices(
        dim_t const nvtxs) noexcept
    {
      setNumRows(nvtxs);
      setNumCols(nvtxs);
    }


    /**
     * @brief Get the number of edges in this graph.
     *
     * @return The number of edges in this graph.
     */
    inline ind_t getNumEdges() const noexcept
    {
      return getNNZ();
    }


    /**
     * @brief Set the number of edges for this graph.
     *
     * @param nedges The number of edges in this graph.
     */
    inline void setNumEdges(
        ind_t const nedges) noexcept
    {
      setNNZ(nedges);
    }


    /**
     * @brief Check if this graph has edge weights present.
     *
     * @return True if this graph has edge weights.
     */
    inline bool hasEdgeWeights() const noexcept
    {
      return m_ewgts;
    }


    /**
     * @brief Set whether this graph has edge weights.
     *
     * @param ewgts Whether this graph has edge weights.
     */
    inline void hasEdgeWeights(
        bool const ewgts) noexcept
    {
      m_ewgts = ewgts;
    }


    /**
     * @brief Get the number of vertex weights in this graph.
     *
     * @return The number of vertex weights in this graph. 
     */
    inline int getNumVertexWeights() const noexcept
    {
      return m_nvwgts;
    }


    /**
     * @brief Set the number of vertex weights in this graph.
     *
     * @param nvwgts The number of vertex weights in this graph.
     */
    inline void setNumVertexWeights(
        int const nvwgts) noexcept
    {
      m_nvwgts = nvwgts;
    }


    /**
     * @brief Go to the beginning of the filestream.
     */
    virtual void resetStream() = 0;


  private:
    /**
     * @brief Indicates whether edge weights are present in this graph file.
     */
    bool m_ewgts;


    /**
     * @brief The number of vertex weights associated with thsi graph file.
     */
    int m_nvwgts;





};




}



#endif
