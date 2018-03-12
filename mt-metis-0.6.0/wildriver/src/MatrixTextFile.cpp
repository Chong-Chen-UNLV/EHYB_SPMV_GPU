/**
 * @file MatrixTextFile.cpp
 * @brief Abstract providing common code for matrix reading/writing classes.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#include "MatrixTextFile.hpp"
#include "MetisFile.hpp"
#include "CSRFile.hpp"




namespace WildRiver
{


/******************************************************************************
* PROTECTED FUNCTIONS *********************************************************
******************************************************************************/


dim_t MatrixTextFile::incRow()
{
  if (current_row >= getNumRows()) {
    throw BadFileStateException(std::string("Attempt to increase current " \
          "row beyond the number of rows: ") + std::to_string(getNumRows()));
  }

  ++current_row;

  return current_row;
}




/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


MatrixTextFile::MatrixTextFile(
    std::string const & name) :
  TextFile(name),
  MatrixFile(),
  read_nnz(0),
  current_row(0)
{
  // do nothing
}


MatrixTextFile::~MatrixTextFile()
{
  // do nothing
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void MatrixTextFile::read(
    ind_t * const rowptr,
    dim_t * const rowind,
    val_t * const rowval,
    double * const progress)
{
  std::vector<matrix_entry_struct> row;

  dim_t nrows = getNumRows();

  dim_t const interval = nrows > 100 ? nrows / 100 : 1;
  double const increment = 1.0/100.0;

  // read in the rows the matrix into our csr
  dim_t j = 0;
  rowptr[0] = j;
  for (dim_t i=0; i<nrows; ++i) {
    dim_t degree;

    dim_t * const rowindStart = rowind+rowptr[i];
    val_t * const rowvalStart = rowval ? rowval+rowptr[i] : nullptr;

    if (!getNextRow(&degree,rowindStart,rowvalStart)) {
      // fewer rows than expected
      throw EOFException(std::string("Failed to read row ") + \
          std::to_string(i) + std::string("/") + std::to_string(nrows)); 
    }

    rowptr[i+1] = rowptr[i]+degree;

    if (progress != nullptr && i % interval == 0) {
      *progress += increment;
    }
  }

  if (rowptr[nrows] != getNNZ()) {
    // we read in the wrong number of non-zeroes
    throw EOFException(std::string("Only found ") + std::to_string(j) + \
        std::string("/") + std::to_string(getNNZ()) + \
        std::string(" non-zeroes in file"));
  }
}


void MatrixTextFile::write(
    ind_t const * const rowptr,
    dim_t const * const rowind,
    val_t const * const rowval)
{
  std::vector<matrix_entry_struct> row;

  dim_t const nrows = getNumRows();

  for (dim_t i=0;i<nrows;++i) {
    // build and insert a new row
    row.clear();
    for (ind_t j=rowptr[i];j<rowptr[i+1];++j) {
      matrix_entry_struct entry;

      entry.ind = rowind[j];
      if (rowval) {
        entry.val = rowval[j];
      } else {
        entry.val = 1;
      }
      row.push_back(entry);
    }
    setNextRow(row);
  }
}




}
