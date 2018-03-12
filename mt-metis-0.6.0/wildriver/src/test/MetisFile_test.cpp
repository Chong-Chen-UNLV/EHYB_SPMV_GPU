/**
 * @file MetisFile_test.cpp
 * @brief Test for reading and writing Metis formatted graphs.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015
 * @version 1
 *
 */




#include <iostream>
#include <fstream>
#include <memory>

#include "MetisFile.hpp"
#include "DomTest.hpp"




using namespace WildRiver;



namespace DomTest
{


static void writeTest(
    std::string const & testFile)
{
  MetisFile graph(testFile);

  graph.setInfo(6,14,0,false);

  wildriver_ind_t xadj[] = {0,2,4,7,10,12,14};
  wildriver_dim_t adjncy[] = {1,2,0,2,0,1,3,2,4,5,3,5,3,4};

  graph.write(xadj,adjncy,nullptr,nullptr);
}


static void readTest(
    std::string const & testFile)
{
  MetisFile graph(testFile);

  wildriver_dim_t nvtxs;
  wildriver_ind_t nedges;
  int nvwgts;
  bool ewgts;

  graph.getInfo(nvtxs,nedges,nvwgts,ewgts);

  testEquals(nvtxs,6);
  testEquals(nedges,14);
  testEquals(nvwgts,0);
  testEquals(ewgts,false);

  std::unique_ptr<wildriver_ind_t[]> xadj(new wildriver_ind_t[nvtxs+1]);
  std::unique_ptr<wildriver_dim_t[]> adjncy(new wildriver_dim_t[nedges]);

  graph.read(xadj.get(),adjncy.get(),nullptr,nullptr,nullptr);

  // test xadj
  testEquals(xadj[0],0);
  testEquals(xadj[1],2);
  testEquals(xadj[2],4);
  testEquals(xadj[3],7);
  testEquals(xadj[4],10);
  testEquals(xadj[5],12);
  testEquals(xadj[6],14);

  // test adjncy
  testEquals(adjncy[0],1);
  testEquals(adjncy[1],2);

  testEquals(adjncy[2],0);
  testEquals(adjncy[3],2);

  testEquals(adjncy[4],0);
  testEquals(adjncy[5],1);
  testEquals(adjncy[6],3);

  testEquals(adjncy[7],2);
  testEquals(adjncy[8],4);
  testEquals(adjncy[9],5);

  testEquals(adjncy[10],3);
  testEquals(adjncy[11],5);

  testEquals(adjncy[12],3);
  testEquals(adjncy[13],4);
}


void Test::run()
{
  std::string testFile("/tmp/metis_test.graph");

  writeTest(testFile);
  readTest(testFile);

}




}
