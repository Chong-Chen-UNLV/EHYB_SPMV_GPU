/**
 * @file MetisFile.hpp
 * @brief Class for reading/writing metis files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_METISFILE_HPP
#define WILDRIVER_METISFILE_HPP




#include "GraphTextFile.hpp"




namespace WildRiver {


class MetisFile : 
  public GraphTextFile
{
  public:
    /**
     * @brief Check if the given filename matches an extension for this 
     * filetype.
     *
     * @param f The filename.
     *
     * @return True if the extension matches this filetype.
     */
    static bool hasExtension(
        std::string const & f);


    /**
     * @brief Close file and free any memory.
     */
    virtual ~MetisFile();


    /**
     * @brief Create a new MetisFile for reading and writing.
     *
     * @param fname The filename/path.
     */
    MetisFile(
        std::string const & fname);


    /**
     * @brief Reset the current position in the matrix file to the first row.
     */
    virtual void firstRow() override;


    /**
     * @brief Reset the current position in the graph to the first vertex.
     */
    virtual void firstVertex() override;


    /**
     * @brief Get the name of this matrix file type.
     *
     * @return The matrix file type name.
     */
    virtual std::string const & getName() const noexcept override
    {
      return NAME;
    } 


    /**
     * @brief Set the next row in the matrix (adjacecny list in the graph).
     *
     * @param numNonZeros The number of non-zeros in the row (output).
     * @param columns The column of each non-zero entry (must be of length at
     * least the number of non-zero entries).
     * @param values The value of each non-zero entry (must be null or of 
     * length at least the number of non-zero entries).
     *
     * @return True if another row was found in the file.
     */
    virtual bool getNextRow(
        dim_t * numNonZeros,
        dim_t * columns,
        val_t * values) override;


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
     * @brief Set the next row in the matrix file.
     *
     * @param next The row to set.
     */
    virtual void setNextRow(
        std::vector<matrix_entry_struct> const & next) override;


    /**
     * @brief Set the adjacency list and vertex weight of the next vertex.
     *
     * @param vwgts The vertex weights for this vertex.
     * @param list The adjacecny list.
     */
    virtual void setNextVertex(
        std::vector<val_t> const & vwgts,
        std::vector<matrix_entry_struct> const & list) override;


  protected:
    /**
     * @brief Determine if teh given line is a comment.
     *
     * @param line The line.
     *
     * @return True if the line is a comment.
     */
    virtual bool isComment(
        std::string const & line) const noexcept override;


    /**
     * @brief Read the header of this matrix file. Populates internal fields
     * with the header information.
     */
    virtual void readHeader() override;


    /**
     * @brief Write the header of this matrix file. The header consists of
     * internal fields set by "setInfo()".
     */
    virtual void writeHeader() override; 


  private:
    /**
     * @brief Name of this filetype.
     */
    static std::string const NAME;


    /**
     * @brief Line buffer.
     */
    std::string m_line;


    /**
     * @brief Get the flags representing the weights associated with this
     * graph.
     *
     * @return 
     */
    int getWeightFlags();




};




}




#endif
