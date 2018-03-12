/**
 * @file GraphTextFile.hpp
 * @brief Abstract class for reading and writing graphs.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */





#ifndef WILDRIVER_GRAPHTEXTFILE_HPP
#define WILDRIVER_GRAPHTEXTFILE_HPP




#include <vector>
#include <string>

#include "GraphFile.hpp"
#include "MatrixTextFile.hpp"




namespace WildRiver
{


class GraphTextFile :
  public TextFile,
  public GraphFile
{
  public:
    /**
     * @brief Create a new graph file.
     *
     * @param fname The filename/path.
     */
    GraphTextFile(
        std::string const & fname);


    /**
     * @brief Destructor.
     */
    virtual ~GraphTextFile();


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
        double * progress) override;


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
        val_t const * adjwgt) override;


    /**
     * @brief Get the filename/path of the graph file.
     *
     * @return The filename/path. 
     */
    virtual std::string const & getFilename() const noexcept
    {
      return TextFile::getFilename();
    }


  protected:
    /**
     * @brief Incremenet the current vertex in the stream.
     *
     * @return The new vertex number. 
     */
    dim_t incVertex();


    /**
     * @brief Go to the beginning of the filestream.
     */
    virtual void resetStream() override 
    {
      TextFile::resetStream();
    }





  private:
    dim_t m_currentVertex;





};




}




#endif
