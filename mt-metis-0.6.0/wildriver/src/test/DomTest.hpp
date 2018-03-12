/**
 * @file DomTest.hpp
 * @brief Top level header for DomTest 
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2015-08-22
 */




#ifndef DOMTEST_HPP
#define DOMTEST_HPP




#include <stdexcept>
#include <string>
#include <iostream>




/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define testEquals(a,b) \
  if ((a) != (b)) { \
    std::cerr << "Test Failed: " << #a << ":'" << a << "' != " << #b << ":'" \
        << b << "' at " << __LINE__ << std::endl; \
    throw TestFailed("Equals Failed"); \
  }


#define testStringEquals(a,b) \
  if ((a).compare(b) != 0) { \
    std::cerr << "Test Failed: " << #a << ":'" << a << "' != " << #b << ":'" \
        << b << "' at " << __LINE__ << std::endl; \
    throw TestFailed("Equals Failed"); \
  }


#define testTrue(a) \
  if (!(a)) { \
    std::cerr << "Test Failed: '" << #a << "' is false at: " << __LINE__ \
      << std::endl; \
    throw TestFailed("True Failed"); \
  }


#define testGreaterThan(a,b) \
  if ((a) <= (b)) { \
    std::cerr << "Test Failed:" << #a << ":'" << a << "' <= " << #b << ":'" \
        << b << "' at " << __LINE__ << std::endl; \
    throw TestFailed("Greater Than Failed"); \
  }


#define testGreaterThanOrEqual(a,b) \
  if ((a) < (b)) { \
    std::cerr << "Test Failed:" << #a << ":'" << a << "' < " << #b << ":'" \
        << b << "' at " << __LINE__ << std::endl; \
    throw TestFailed("Greater Than Or Equal Failed"); \
  }






namespace DomTest
{


class TestFailed : public std::logic_error
{
  public:
    TestFailed(std::string const & str) :
        std::logic_error(str)
    {
    }
};


/******************************************************************************
* TEST CLASS ******************************************************************
******************************************************************************/


static void removeFile(
    std::string const & file)
{
  remove(file.c_str());
}


class Test
{
  protected:
    void run();


  public:
    Test()
    {
    };

    bool evaluate()
    {
      try {
        run();
      } catch (std::exception & e) {
        std::cerr << "FAILED: " << e.what() << std::endl;
        // fail
        return false;
      }

      // pass
      return true;
    };

    virtual ~Test() 
    {
    }
};




}




/******************************************************************************
* MAIN ************************************************************************
******************************************************************************/


int main(
    int argc,
    char ** argv)
{
  // make sure we don't have any useless arguments
  if (argc > 1) {
    for (int i=1; i<argc; ++i) {
      std::cerr << "Unused parameter: " << argv[i] << std::endl;
    }
    return 0;
  }

  DomTest::Test test;

  if (test.evaluate()) {
    return 0;
  } else {
    return 1;
  }
}




#endif
