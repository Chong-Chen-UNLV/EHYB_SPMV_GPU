/**
 * @file graph_test.c
 * @brief Unit tests for the graph.h file.
 * @author Dominique LaSalle <dominique@domnet.org>
 * Copyright 2016
 * @version 1
 * @date 2016-07-31
 */




#include "test.h"
#include "graph.h"




int test(void)
{
  /* test graph_create */
  graph_type * graph = graph_create(4);
  TESTEQUALS(graph->level,(size_t)0,"%zu");
  TESTEQUALS(graph->nvtxs,0,"%"PF_VTX_T);
  TESTEQUALS(graph->nedges,0,"%"PF_ADJ_T);
  TESTEQUALS(graph->mincut,0,"%"PF_WGT_T);
  TESTEQUALS(graph->minvol,0,"%"PF_VTX_T);
  TESTEQUALS(graph->ngroup,0,"%"PF_PID_T);
  TESTEQUALS((void*)graph->group,NULL,"%p");
  TESTTRUE(graph->mynvtxs != NULL);
  TESTTRUE(graph->mynedges != NULL);
  TESTTRUE(graph->xadj != NULL);
  TESTTRUE(graph->vwgt != NULL);
  TESTTRUE(graph->adjncy != NULL);
  TESTTRUE(graph->adjwgt != NULL);
  TESTEQUALS((void*)graph->label,NULL,"%p");
  TESTEQUALS((void*)graph->rename,NULL,"%p");
  TESTEQUALS((void*)graph->cmap,NULL,"%p");
  TESTEQUALS((void*)graph->nislands,NULL,"%p");
  TESTEQUALS(graph->uniformvwgt,0,"%"PF_WGT_T);
  TESTEQUALS(graph->uniformadjwgt,0,"%"PF_WGT_T);
  TESTEQUALS(graph->tvwgt,(twgt_type)0,"%"PF_TWGT_T);
  TESTEQUALS(graph->invtvwgt,0.0,"%"PF_REAL_T);
  TESTEQUALS((void*)graph->where,NULL,"%p");
  TESTEQUALS((void*)graph->pwgts,NULL,"%p");
  TESTEQUALS((void*)graph->vsinfo,NULL,"%p");
  TESTEQUALS((void*)graph->esinfo,NULL,"%p");
  TESTEQUALS((void*)graph->kwinfo,NULL,"%p");
  TESTEQUALS(graph->free_xadj,1,"%d");
  TESTEQUALS(graph->free_vwgt,1,"%d");
  TESTEQUALS(graph->free_adjncy,1,"%d");
  TESTEQUALS(graph->free_adjwgt,1,"%d");
  TESTEQUALS((void*)graph->coarser,NULL,"%p");
  TESTEQUALS((void*)graph->finer,NULL,"%p");
  

  /* test graph_setup */

  /* test graph_distribute */

  /* test graph_gather */

  /* test graph_setup_coarse */

  /* test graph_setup_twgts */

  /* test graph_alloc_partmemory */

  /* test graph_free */

  /* test graph_free_rdata */

  /* test graph_imbalance */

  /* test graph_imbalance_diff */

  /* test graph_cut */

  /* test graph_isbalanced */

  /* test graph_readjust_memory */

  /* test graph_extract_halves */

  /* test graph_size */

  /* test graph_calc_dist */

  return 0;
}
