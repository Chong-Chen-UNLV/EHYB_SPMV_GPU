/**
 * @file MatrixGraphFile.hpp
 * @brief An adapter class for treating matrices as graphs.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_MATRIXGRAPHFILE_HPP
#define WILDRIVER_MATRIXGRAPHFILE_HPP




#include <vector>
#include <string>
#include <memory>

#include "IGraphFile.hpp"
#include "IMatrixFile.hpp"




namespace WildRiver
{


class MatrixGraphFile :
  public IGraphFile
{
  public:
    /**
     * @brief Create an adapter so that a matrix can be treated as a graph. The
     * supplied pointer will be free'd when this object is destructed.
     *
     * @param mfile A pointer the MatrixFile to adapt.
     */
    MatrixGraphFile(
        std::unique_ptr<IMatrixFile>& file);


    /**
     * @brief Frees the underlying matrix file.
     */
    ~MatrixGraphFile();


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
    void setInfo(
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
    void read(
        ind_t * xadj,
        dim_t * adjncy,
        val_t * vwgt,
        val_t * adjwgt,
        double * progress) override;


    /**
     * @brief Write a graph file from the given CSR structure.
     *
     * @param xadj The adjacency list pointer (length nvtxs+1).
     * @param adjncy The adjacency list (length nedges).
     * @param vwgt The vertex weights. 
     * @param adjwgt The edge weights.
     */
    void write(
        ind_t const * xadj,
        dim_t const * adjncy,
        val_t const * vwgt,
        val_t const * adjwgt) override;


    /**
     * @brief Set the cursor to the first vertex in the file.
     */
    void firstVertex() override;


    /**
     * @brief Get the information of the next vertex.
     *
     * @param vertexWeights The vertex weight(s) (must be null or at least of
     * length equal to the number of constraints).
     * @param numEdges The number of edges incident to this vertex (output).
     * @param edgeDests The destination of each edge leaving this vertex (must
     * be of length equal to the number of edges of the vertex).
     * @param edgeWeights The weight of each edge leaving this vertex (must
     * null or be of length equal to the number of edges of the vertex).
     *
     * @return True if another vertex was found in the file.
     */
    virtual bool getNextVertex(
        val_t * vertexWeights,
        dim_t * numEdges,
        dim_t * edgeDests,
        val_t * edgeWeights) override;


    /**
     * @brief Set the adjacency list and vertex weight of the next vertex.
     *
     * @param vwgts The vertex weights for this vertex.
     * @param list The adjacecny list.
     */
    void setNextVertex(
        std::vector<val_t> const & vwgts,
        std::vector<matrix_entry_struct> const & list) override;


    /**
     * @brief Get the name of the filetype.
     *
     * @return The name of the filetype.
     */
    std::string const & getName() const noexcept override;


    /**
     * @brief Get the filename/path of the graph/matrix file.
     *
     * @return The filename/path. 
     */
    virtual std::string const & getFilename() const noexcept override
    {
      return m_file->getFilename();
    }


  private:
    std::unique_ptr<IMatrixFile> m_file;


    // disable copying
    MatrixGraphFile(
        MatrixGraphFile const & rhs);
    MatrixGraphFile& operator=(
        MatrixGraphFile const & rhs);



};


}




#endif
