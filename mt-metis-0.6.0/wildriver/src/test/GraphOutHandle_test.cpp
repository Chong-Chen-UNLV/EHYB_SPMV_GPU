/**
 * @file GraphOutHandle_test.cpp
 * @brief Test for reading matrices
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015
 * @version 1
 *
 */




#include <iostream>
#include <fstream>
#include <memory>

#include "GraphOutHandle.hpp"
#include "DomTest.hpp"




static size_t const NVTXS = 6;
static size_t const NEDGES = 7;




using namespace WildRiver;




namespace DomTest
{


static void writeGraph(
    std::string const & testFile)
{
  wildriver_dim_t nvtxs = 6;
  wildriver_ind_t nedges = 14;
  wildriver_ind_t xadj[] = {0,2,4,7,10,12,14};
  wildriver_dim_t adjncy[] = {1,2,0,2,0,1,3,2,4,5,3,5,3,4};
  wildriver_val_t adjwgt[] = {1,2,3,4,5,6,7,8,9,1,2,3,4,5};


  GraphOutHandle handle(testFile);

  handle.setInfo(nvtxs,nedges,0,true);

  handle.writeGraph(xadj,adjncy,NULL,adjwgt);
}


static void readGraph(
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
  std::string metisFile("/tmp/GraphOutHandle_test.graph");
  writeGraph(metisFile);
  readMetis(metisFile);

  // generate test csr file
  std::string csrFile("/tmp/GraphOutHandle_test.csr");
  writeGraph(csrFile);
  readGraph(csrFile);
}




}
