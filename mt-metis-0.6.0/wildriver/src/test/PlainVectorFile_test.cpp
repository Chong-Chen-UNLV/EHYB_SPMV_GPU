/**
 * @file PlainVectorFile_test.cpp
 * @brief 
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015
 * @version 1
 * @date 2016-02-14
 */




#include <cstdio>
#include <iostream>
#include <fstream>
#include <memory>

#include "PlainVectorFile.hpp"
#include "DomTest.hpp"




using namespace WildRiver;




namespace DomTest
{


static void writeTest(
    std::string const & testFile)
{
  removeFile(testFile);

  PlainVectorFile vec(testFile);

  vec.setSize(10);

  testEquals(vec.getSize(),10);

  wildriver_val_t vals[] = {1,2,3,4,5,6,7,8,9,0};

  vec.write(vals,nullptr);
}


static void readTest(
    std::string const & testFile)
{
  PlainVectorFile vec(testFile);

  dim_t n = vec.getSize();
  testEquals(n,10);

  std::unique_ptr<wildriver_val_t[]> vals(new wildriver_val_t[n]);

  vec.read(vals.get(),nullptr);

  testEquals(vals[0],1);
  testEquals(vals[1],2);
  testEquals(vals[2],3);
  testEquals(vals[3],4);
  testEquals(vals[4],5);
  testEquals(vals[5],6);
  testEquals(vals[6],7);
  testEquals(vals[7],8);
  testEquals(vals[8],9);
  testEquals(vals[9],0);
}


void Test::run()
{
  std::string testFile("/tmp/test.vec");

  writeTest(testFile);
  readTest(testFile);

  removeFile(testFile);
}




}
