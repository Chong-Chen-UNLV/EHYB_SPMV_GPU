/**
 * @file MetisFile.cpp
 * @brief Class functions for reading and writing metis graph files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#include "MetisFile.hpp"




namespace WildRiver 
{


/******************************************************************************
* INTERNAL TYPES **************************************************************
******************************************************************************/


enum {
  HAS_NOWEIGHTS = 0,
  HAS_EDGEWEIGHTS = 1,
  HAS_VERTEXWEIGHTS = 10
};




/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/


namespace
{


static size_t const BUFFER_SIZE = 4096;


}


std::string const MetisFile::NAME = "Metis";




/******************************************************************************
* PUBLIC STATIC FUNCTIONS *****************************************************
******************************************************************************/


bool MetisFile::hasExtension(
    std::string const & f)
{
  std::vector<std::string> extensions;

  extensions.push_back(".graph");
  extensions.push_back(".metis");
  extensions.push_back(".chaco");

  return TextFile::matchExtension(f,extensions);
}




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


int MetisFile::getWeightFlags()
{
  int flags = HAS_NOWEIGHTS;

  if (hasEdgeWeights()) {
    flags |= HAS_EDGEWEIGHTS;
  }

  if (getNumVertexWeights() > 0) {
    flags |= HAS_VERTEXWEIGHTS;
  }

  return flags;
}




/******************************************************************************
* PROTECTED FUNCTIONS *********************************************************
******************************************************************************/


bool MetisFile::isComment(
    std::string const & line) const noexcept
{
  // awful solution since I can't declare this statically in c++ -- at
  // somepoint generate all 256 entries using template programming
  bool comment_chars[256] = {false};
  comment_chars['#'] = true;
  comment_chars['%'] = true;
  comment_chars['"'] = true;
  comment_chars['/'] = true;

  return line.size() > 0 && comment_chars[static_cast<uint8_t>(line[0])];
}


void MetisFile::readHeader()
{
  if (!isOpenRead()) {
    // open our file for reading if not already
    openRead();
  }

  // get the first line
  std::string line;
  nextNoncommentLine(line);

  // parse out my header
  size_t offset;
  dim_t nvtxs = std::stoull(line,&offset,10);
  line = line.substr(offset);

  setNumVertices(nvtxs);

  ind_t nedges = std::stoull(line,&offset,10)*2;
  line = line.substr(offset);

  setNumEdges(nedges);

  // handle weights
  if (line.find_first_not_of(" \t") != std::string::npos) {
    int flags = std::stoul(line,&offset,10);  
    line = line.substr(offset);

    if (flags & HAS_EDGEWEIGHTS) {
      hasEdgeWeights(true);
    }

    if (flags & HAS_VERTEXWEIGHTS) {
      setNumVertexWeights(std::stoul(line,&offset,10));
    }
  }
}


void MetisFile::writeHeader()
{
  dim_t nvtxs = getNumVertices();
  ind_t nedges = getNumEdges();

  // make sure we have a valid number of edges/nnz
  if (nedges % 2 != 0) {
    throw BadParameterException("Metis files are required to be symmetric: " \
        "odd number of non-zeroes specified.");
  }

  // set header to write mode if not already
  if (!isOpenWrite()) {
    openWrite();
  }

  // write the header -- edges in metis files are undirected (symmetric).
  getStream() << nvtxs << " " << (nedges/2);

  int weightflag = getWeightFlags();  

  if (weightflag != HAS_NOWEIGHTS) {
    // write weight flags
    getStream() << " " << weightflag;
    if (weightflag & HAS_VERTEXWEIGHTS) {
      // write number of vertex weights
      getStream() << " " << getNumVertexWeights();
    }
  } 

  getStream() << std::endl;
}




/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS **************************************************
******************************************************************************/


MetisFile::MetisFile(
    std::string const & fname) : 
  GraphTextFile(fname),
  m_line(BUFFER_SIZE,'\0')
{
}


MetisFile::~MetisFile()
{
}




/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


void MetisFile::firstRow()
{
  resetStream();

  // go past all comments -- read and discard header
  std::string line;
  nextNoncommentLine(line);
}


void MetisFile::firstVertex()
{
  firstRow();
}


bool MetisFile::getNextRow(
    dim_t * const numNonZeros,
    dim_t * const columns,
    val_t * const values)

{
  return getNextVertex(nullptr,numNonZeros,columns,values);
}


bool MetisFile::getNextVertex(
        val_t * vertexWeights,
        dim_t * numEdges,
        dim_t * edgeDests,
        val_t * edgeWeights)
{
  dim_t const ncon = getNumVertexWeights();

  // get my m_line
  if (!nextNoncommentLine(m_line)) {
    return false;
  }

  char * sptr;
  char * eptr = (char*)m_line.data();

  // read in vertex weights
  for (dim_t k=0; k<ncon; ++k) {
    sptr = eptr;
    val_t val = std::strtod(sptr,&eptr);
    if (sptr == eptr) {
      throw BadFileException(std::string("Failed to read vertex weight on " \
            "line ") + std::to_string(getCurrentLine()));
    }
    if (vertexWeights) {
      vertexWeights[k] = val;
    }
  }

  dim_t const nvtxs = getNumVertices();
  bool ewgts = hasEdgeWeights();

  dim_t degree = 0;
  // read in edges
  while (true) {
    dim_t dst;
    val_t wgt;

    sptr = eptr;
    dst = static_cast<dim_t>(std::strtoull(sptr,&eptr,10))-1;
    if (sptr == eptr) {
      break;
    }

    // make sure this is a valid edge
    if (dst >= nvtxs) {
      throw BadFileException(std::string("Edge with destination of ") + \
          std::to_string(dst) + std::string("/") + \
          std::to_string(nvtxs));
    }

    if (ewgts) {
      sptr = eptr;
      wgt = static_cast<val_t>(std::strtod(sptr,&eptr));
      if (sptr == eptr) {
        throw BadFileException(std::string("Could not read edge weight at "
              "line ") + std::to_string(getCurrentLine()));
      }
    }

    edgeDests[degree] = dst;
    if (edgeWeights) {
      if (ewgts) {
        edgeWeights[degree] = wgt;
      } else {
        edgeWeights[degree] = static_cast<val_t>(1);
      }
    }
    
    ++degree;
  }

  *numEdges = degree;

  // indicate that we successfully found a vertex
  return true;
}


void MetisFile::setNextRow(
    std::vector<matrix_entry_struct> const & row)
{
  dim_t const nadj = row.size();
  bool const ewgts = hasEdgeWeights();

  for (dim_t j=0; j<nadj; ++j) {
    matrix_entry_struct e = row[j];
    getStream() << (e.ind+1);
    if (ewgts) {
      getStream() << " " << e.val;
    }
    if (j < nadj-1) {
      // add space
      getStream() << " ";
    }
  }
  getStream() << std::endl;

  incVertex();
}


void MetisFile::setNextVertex(
    std::vector<val_t> const & vwgts,
    std::vector<matrix_entry_struct> const & list)
{
  dim_t const ncon = getNumVertexWeights();
  dim_t const nadj = list.size();

  // set vertex weights
  for (dim_t k=0; k<ncon; ++k) {
    getStream() << vwgts[k];
    if (k < ncon-1 || nadj > 0) {
      getStream() << " ";
    }
  }

  // set adjacency list and edge weights
  setNextRow(list);
}




}
