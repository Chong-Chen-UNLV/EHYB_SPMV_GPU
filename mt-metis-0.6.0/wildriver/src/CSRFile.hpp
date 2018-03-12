/**
 * @file CSRFile.hpp
 * @brief Class for reading/writing metis files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_CSRFILE_HPP
#define WILDRIVER_CSRFILE_HPP




#include "MatrixTextFile.hpp"




namespace WildRiver {


class CSRFile : 
  public MatrixTextFile
{
  public:
    /**
     * @brief Name of this filetype.
     */
    static std::string const NAME;


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
     * @brief Create a new CSRFile for reading and writing.
     *
     * @param fname The filename/path.
     */
    CSRFile(
        std::string const & fname);


    /**
     * @brief Close file and free any memory.
     */
    virtual ~CSRFile();


    /**
     * @brief Reset the current position in the matrix file to the first row.
     */
    virtual void firstRow() override;


    /**
     * @brief Get the next row in the matrix (adjacecny list in the graph).
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
     * @brief Set the next row in the matrix file.
     *
     * @param next The row to set.
     */
    virtual void setNextRow(
        std::vector<matrix_entry_struct> const & next) override;


    /**
     * @brief Get the name of this matrix file type.
     *
     * @return The matrix file type name.
     */
    virtual std::string const & getName() const noexcept override
    {
      return NAME;
    } 


  protected:
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


    /**
     * @brief Determine the given line is a comment.
     *
     * @param line The line.
     *
     * @return True if the line is a comment.
     */
    virtual bool isComment(
        std::string const & line) const noexcept override;


  private:
    std::string m_line;




};




}




#endif
