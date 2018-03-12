/**
 * @file PlainVectorFile.cpp
 * @brief Implemenation of PlainVectorFile class.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-08
 */





#include "PlainVectorFile.hpp"




namespace WildRiver
{


/******************************************************************************
* PUBLIC STATIC FUNCTIONS *****************************************************
******************************************************************************/


bool PlainVectorFile::hasExtension(
    std::string const & f)
{
  std::vector<std::string> extensions;

  extensions.push_back(".txt");
  extensions.push_back(".vec");
  extensions.push_back(".perm");
  extensions.push_back(".part");
  extensions.push_back(".cluster");

  return TextFile::matchExtension(f,extensions);
}




/******************************************************************************
* PROTECTED FUNCTIONS *********************************************************
******************************************************************************/


bool PlainVectorFile::isComment(
    std::string const & line) const noexcept
{
  if (line.size() > 0) {
    switch (line[0]) {
      case '#':
      case '%':
      case '/':
        return true;
      default:
        return false;
    }
  } else {
    return false;
  }
}


/******************************************************************************
* CONSTRUCTORS / DESTRUCTOR ***************************************************
******************************************************************************/


PlainVectorFile::PlainVectorFile(
    std::string const & name) :
  VectorTextFile(name),
  m_buffer()
{
  // do nothing
}


PlainVectorFile::~PlainVectorFile()
{
  // do nothing
}


/******************************************************************************
* PUBLIC FUNCTIONS ************************************************************
******************************************************************************/


ind_t PlainVectorFile::getSize()
{
  if (!isSizeSet()) {
    size_t nlines = 0;
    std::string m_buffer;

    if (!isOpenRead()) {
      openRead();
    }

    // cout non-comment lines
    while (nextNoncommentLine(m_buffer)) {
      ++nlines;
    }

    VectorFile::setSize(nlines);

    // go back to the beginning
    resetStream();
  }

  return VectorFile::getSize();
}


void PlainVectorFile::read(
    val_t * const vals,
    double * progress)
{
  if (!isOpenRead()) {
    openRead();
  }

  const size_t n = getSize();

  resetStream();

  ind_t const interval = n > 100 ? n / 100 : 1;
  double const increment = 1.0 / 100.0;
  
  for (size_t i = 0; i < n; ++i) {
    if (!nextNoncommentLine(m_buffer)) {
      throw EOFException("Hit end of file before getting next value");
    }

    if (typeid(val_t) == typeid(double) || typeid(val_t) == typeid(float)) {
      vals[i] = static_cast<val_t>(std::stod(m_buffer));
    } else {
      vals[i] = static_cast<val_t>(std::stoll(m_buffer));
    }

    if (progress != nullptr && i % interval == 0) {
      *progress += increment;
    }
  }
}


void PlainVectorFile::write(
    val_t const * const vals,
    double  * progress)
{
  if (!isSizeSet()) {
    throw UnsetInfoException("Size of vector is not set before call to " \
        "write()");
  }

  if (!isOpenWrite()) {
    openWrite();
  }

  const size_t n = getSize();
  for (size_t i = 0; i < n; ++i) {
	  getStream() << vals[i] << std::endl;	
  }
}




}
