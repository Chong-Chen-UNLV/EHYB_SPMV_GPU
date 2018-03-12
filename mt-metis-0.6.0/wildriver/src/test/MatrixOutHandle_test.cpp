/**
 * @file MatrixOutHandle_test.cpp
 * @brief Test for reading matrices
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015
 * @version 1
 *
 */




#include <iostream>
#include <fstream>
#include <memory>

#include "MatrixOutHandle.hpp"
#include "DomTest.hpp"




static size_t const NVTXS = 6;
static size_t const NEDGES = 7;




using namespace WildRiver;




namespace DomTest
{


static void writeSparse(
    std::string const & testFile)
{
  wildriver_dim_t nrows = 6, ncols = 6;
  wildriver_ind_t nnz = 14;
  wildriver_ind_t rowptr[] = {0,2,4,7,10,12,14};
  wildriver_dim_t rowind[] = {1,2,0,2,0,1,3,2,4,5,3,5,3,4};
  wildriver_val_t rowval[] = {1,2,3,4,5,6,7,8,9,1,2,3,4,5};


  MatrixOutHandle handle(testFile);

  handle.setInfo(nrows,ncols,nnz);

  handle.writeSparse(rowptr,rowind,rowval);
}


static void readSparse(
    std::string const & testFile)
{
  std::fstream stream(testFile,std::fstream::in);
  std::string line;

  // test matrix itself
  std::getline(stream,line);
  testStringEquals(line,"2 1 3 2"); 

  std::getline(stream,line);
  testStringEquals(line,"1 3 3 4"); 

  std::getline(stream,line);
  testStringEquals(line,"1 5 2 6 4 7");

  std::getline(stream,line);
  testStringEquals(line,"3 8 5 9 6 1");

  std::getline(stream,line);
  testStringEquals(line,"4 2 6 3");

  std::getline(stream,line);
  testStringEquals(line,"4 4 5 5");
}


static void readMetis(
    std::string const & testFile)
{
  std::fstream stream(testFile,std::fstream::in);
  std::string line;

  // test header
  std::getline(stream,line);
  testStringEquals(line,"6 7 1");
  

  // test the graph itself
  std::getline(stream,line);
  testStringEquals(line,"2 1 3 2"); 

  std::getline(stream,line);
  testStringEquals(line,"1 3 3 4"); 

  std::getline(stream,line);
  testStringEquals(line,"1 5 2 6 4 7");

  std::getline(stream,line);
  testStringEquals(line,"3 8 5 9 6 1");

  std::getline(stream,line);
  testStringEquals(line,"4 2 6 3");

  std::getline(stream,line);
  testStringEquals(line,"4 4 5 5");
}


void Test::run()
{
  // generate test metis file
  std::string metisFile("/tmp/MatrixOutHandle_test.graph");
  writeSparse(metisFile);
  readMetis(metisFile);

  // generate test csr file
  std::string csrFile("/tmp/MatrixOutHandle_test.csr");
  writeSparse(csrFile);
  readSparse(csrFile);
}




}
