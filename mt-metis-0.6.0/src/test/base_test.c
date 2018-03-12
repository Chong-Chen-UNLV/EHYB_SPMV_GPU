/**
 * @file base_test.c
 * @brief Unit tests for the base.h file.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2016
 * @version 1
 * @date 2016-07-10
 */




#include "test.h"
#include "graph.h"
#include "base.h"




int test(void)
{
  graph_type graph;
  graphdist_type dist;
  dist.nthreads = 7;
  dist.shift = 5;
  dist.offset = 1 << dist.shift;
  dist.mask = dist.offset - 1;

  graph.dist = dist;

  /* gvtx_to_lvtx test */
  TESTEQUALS(0,gvtx_to_lvtx(dist.offset,dist),"%"PF_VTX_T);
  TESTEQUALS(0,gvtx_to_lvtx(dist.offset*2,dist),"%"PF_VTX_T);
  TESTEQUALS(5,gvtx_to_lvtx(dist.offset*2+5,dist),"%"PF_VTX_T);

  /* lvtx_to_gvtx test */
  TESTEQUALS(dist.offset,lvtx_to_gvtx(0,0,dist),"%"PF_VTX_T);
  TESTEQUALS(2*dist.offset,lvtx_to_gvtx(0,1,dist),"%"PF_VTX_T);
  TESTEQUALS(dist.offset+10,lvtx_to_gvtx(10,0,dist),"%"PF_VTX_T);

  /* gvtx_to_tid test */
  TESTEQUALS(0,gvtx_to_tid(dist.offset,dist),"%"PF_TID_T);
  TESTEQUALS(1,gvtx_to_tid(dist.offset*2,dist),"%"PF_TID_T);
  TESTEQUALS(2,gvtx_to_tid(dist.offset*3+10,dist),"%"PF_TID_T);

  /* max_gvtx test */
  TESTEQUALS((dist.nthreads+1)*(dist.offset),max_gvtx(&graph),"%"PF_VTX_T);

  return 0; 
}


