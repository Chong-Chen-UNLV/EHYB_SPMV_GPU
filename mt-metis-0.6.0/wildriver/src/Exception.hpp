/**
 * @file Exception.hpp
 * @brief Exceptions used within WildRiver
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_EXCEPTION_HPP
#define WILDRIVER_EXCEPTION_HPP




#include <stdexcept>




namespace WildRiver
{


class BadParameterException : public std::runtime_error 
{
  public:
    BadParameterException(
        std::string const & str) : 
      std::runtime_error(str)
    {
    }
};


class UnknownExtensionException : public std::runtime_error 
{
  public:
    UnknownExtensionException(
        std::string const & str) : 
      std::runtime_error(str)
    {
    }
};


class BadFileException : public std::runtime_error 
{
  public:
    BadFileException(
        std::string const & str) : 
      std::runtime_error(str)
    {
    }
};


class EOFException : public BadFileException 
{
  public:
    EOFException(
        std::string const & str) : 
      BadFileException(str)
    {
    }
};


class BadFileStateException : public std::logic_error
{
  public:
    BadFileStateException(
        std::string const & str) :
      std::logic_error(str)
  {
  }
};


class UnsetInfoException : public BadFileStateException
{
  public:
    UnsetInfoException(
        std::string const & str) :
      BadFileStateException(str)
  {
  }
};


class OutOfMemoryException : public std::runtime_error 
{
  public:
    OutOfMemoryException(
        size_t const nbytes) : 
      std::runtime_error("Allocation failed: " + \
          std::to_string(nbytes) + " bytes")
    {
    }
};





}


#endif
