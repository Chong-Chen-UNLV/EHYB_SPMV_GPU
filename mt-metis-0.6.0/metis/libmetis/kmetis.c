/*!
\file  
\brief The top-level routines for  multilevel k-way partitioning that minimizes
       the edge cut.

\date   Started 7/28/1997
\author George  
\author Copyright 1997-2011, Regents of the University of Minnesota 
\version\verbatim $Id: kmetis.c 17715 2014-10-02 16:08:42Z dominique $ \endverbatim
*/

#include "metislib.h"


/*************************************************************************/
/*! This function is the entry point for MCKMETIS */
/*************************************************************************/
int METIS_PartGraphKway(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy, 
          idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *nparts, 
          real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *objval, 
          idx_t *part)
{
  int sigrval=0, renumber=0;
  graph_t *graph;
  ctrl_t *ctrl;

  /* set up malloc cleaning code and signal catchers */
  if (!gk_malloc_init()) 
    return METIS_ERROR_MEMORY;

  gk_sigtrap();

  if ((sigrval = gk_sigcatch()) != 0)
    goto SIGTHROW;


  /* set up the run parameters */
  ctrl = SetupCtrl(METIS_OP_KMETIS, options, *ncon, *nparts, tpwgts, ubvec);
  if (!ctrl) {
    gk_siguntrap();
    return METIS_ERROR_INPUT;
  }

  /* if required, change the numbering to 0 */
  if (ctrl->numflag == 1) {
    Change2CNumbering(*nvtxs, xadj, adjncy);
    renumber = 1;
  }

  /* set up the graph */
  graph = SetupGraph(ctrl, *nvtxs, *ncon, xadj, adjncy, vwgt, vsize, adjwgt);

  /* set up multipliers for making balance computations easier */
  SetupKWayBalMultipliers(ctrl, graph);

  /* set various run parameters that depend on the graph */
  if (ctrl->iptype == METIS_IPTYPE_METISRB) {
    ctrl->CoarsenTo = gk_max((*nvtxs)/(40*gk_log2(*nparts)), 30*(*nparts));
    ctrl->CoarsenTo = 10*(*nparts);
    ctrl->nIparts   = (ctrl->CoarsenTo == 30*(*nparts) ? 4 : 5);
  }
  else {
    ctrl->CoarsenTo = 10*(*nparts);
    ctrl->nIparts   = 10;
  }

  /* take care contiguity requests for disconnected graphs */
  if (ctrl->contig && !IsConnected(graph, 0)) 
    gk_errexit(SIGERR, "METIS Error: A contiguous partition is requested for a non-contiguous input graph.\n");
    
  /* allocate workspace memory */  
  AllocateWorkSpace(ctrl, graph);

  /* start the partitioning */
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startwctimer(ctrl->TotalTmr));

  *objval = MlevelKWayPartitioning(ctrl, graph, part);

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopwctimer(ctrl->TotalTmr));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));

  /* clean up */
  FreeCtrl(&ctrl);

SIGTHROW:
  /* if required, change the numbering back to 1 */
  if (renumber)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);

  gk_siguntrap();
  gk_malloc_cleanup(0);

  return metis_rcode(sigrval);
}


/*************************************************************************/
/*! This function computes a k-way partitioning of a graph that minimizes
    the specified objective function.

    \param ctrl is the control structure
    \param graph is the graph to be partitioned
    \param part is the vector that on return will store the partitioning

    \returns the objective value of the partitoning. The partitioning 
             itself is stored in the part vector.
*/
/*************************************************************************/
idx_t MlevelKWayPartitioning(ctrl_t *ctrl, graph_t *graph, idx_t *part)
{
  idx_t i, objval=0, curobj=0, bestobj=0;
  real_t curbal=0.0, bestbal=0.0;
  graph_t *cgraph;


  for (i=0; i<ctrl->ncuts; i++) {
    cgraph = CoarsenGraph(ctrl, graph);

    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startwctimer(ctrl->InitPartTmr));
    AllocateKWayPartitionMemory(ctrl, cgraph);

    /* Compute the initial partitioning */
    switch (ctrl->iptype) {
      case METIS_IPTYPE_METISRB:
        FreeWorkSpace(ctrl); /* Release the work space, for the recursive metis call */
        InitKWayPartitioningRB(ctrl, cgraph);
        AllocateWorkSpace(ctrl, graph); /* Re-allocate the work space */
        break;

      case METIS_IPTYPE_GROW:
        AllocateRefinementWorkSpace(ctrl, 2*cgraph->nedges);
        InitKWayPartitioningGrow(ctrl, cgraph);
        break;

      default:
        gk_errexit(SIGERR, "Unknown iptype: %d\n", ctrl->iptype);
    }

    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopwctimer(ctrl->InitPartTmr));
    IFSET(ctrl->dbglvl, METIS_DBG_IPART, printf("Initial %"PRIDX \
          "-way partitioning cut: %"PRIDX"\n", ctrl->nparts, objval));

    RefineKWay(ctrl, graph, cgraph);

    switch (ctrl->objtype) {
      case METIS_OBJTYPE_CUT:
        curobj = graph->mincut;
        break;

      case METIS_OBJTYPE_VOL:
        curobj = graph->minvol;
        break;

      default:
        gk_errexit(SIGERR, "Unknown objtype: %d\n", ctrl->objtype);
    }

    curbal = ComputeLoadImbalanceDiff(graph, ctrl->nparts, ctrl->pijbm, ctrl->ubfactors);

    if (i == 0 
        || (curbal <= 0.0005 && bestobj > curobj)
        || (bestbal > 0.0005 && curbal < bestbal)) {
      icopy(graph->nvtxs, graph->where, part);
      bestobj = curobj;
      bestbal = curbal;
    }

    FreeRData(graph);

    if (bestobj == 0)
      break;
  }

  FreeGraph(&graph);

  return bestobj;
}


/*************************************************************************/
/*! This function computes the initial k-way partitioning using PMETIS 
*/
/*************************************************************************/
void InitKWayPartitioning(ctrl_t *ctrl, graph_t *graph)
{
  switch (ctrl->iptype) {
    case METIS_IPTYPE_METISRB:
      InitKWayPartitioningRB(ctrl, graph);
      break;

    case METIS_IPTYPE_GROW:
      InitKWayPartitioningGrow(ctrl, graph);
      break;

    default:
      gk_errexit(SIGERR, "Unknown iptype: %d\n", ctrl->iptype);
  }
}


/*************************************************************************/
/*! Compute the initial k-way partitioning using PMETIS */
/*************************************************************************/
void InitKWayPartitioningRB(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, options[METIS_NOPTIONS], curobj=0;
  idx_t *bestwhere=NULL;
  real_t *ubvec=NULL;
  int status;

  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NITER]   = 10;
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_NO2HOP]  = ctrl->no2hop;
  options[METIS_OPTION_ONDISK]  = ctrl->ondisk;


  ubvec = rmalloc(graph->ncon, "InitKWayPartitioning: ubvec");
  for (i=0; i<graph->ncon; i++) 
    ubvec[i] = (real_t)pow(ctrl->ubfactors[i], 1.0/log(ctrl->nparts));


  switch (ctrl->objtype) {
    case METIS_OBJTYPE_CUT:
    case METIS_OBJTYPE_VOL:
      options[METIS_OPTION_NCUTS] = ctrl->nIparts;
      status = METIS_PartGraphRecursive(&graph->nvtxs, &graph->ncon, 
                   graph->xadj, graph->adjncy, graph->vwgt, graph->vsize, 
                   graph->adjwgt, &ctrl->nparts, ctrl->tpwgts, ubvec, 
                   options, &curobj, graph->where);

      if (status != METIS_OK)
        gk_errexit(SIGERR, "Failed during initial partitioning\n");

      break;

#ifdef XXX /* This does not seem to help */
    case METIS_OBJTYPE_VOL:
      bestwhere = imalloc(graph->nvtxs, "InitKWayPartitioning: bestwhere");
      options[METIS_OPTION_NCUTS] = 2;

      ntrials = (ctrl->nIparts+1)/2;
      for (i=0; i<ntrials; i++) {
        status = METIS_PartGraphRecursive(&graph->nvtxs, &graph->ncon, 
                     graph->xadj, graph->adjncy, graph->vwgt, graph->vsize, 
                     graph->adjwgt, &ctrl->nparts, ctrl->tpwgts, ubvec, 
                     options, &curobj, graph->where);
        if (status != METIS_OK)
          gk_errexit(SIGERR, "Failed during initial partitioning\n");

        curobj = ComputeVolume(graph, graph->where);

        if (i == 0 || bestobj > curobj) {
          bestobj = curobj;
          if (i < ntrials-1)
            icopy(graph->nvtxs, graph->where, bestwhere);
        }

        if (bestobj == 0)
          break;
      }
      if (bestobj != curobj)
        icopy(graph->nvtxs, bestwhere, graph->where);

      break;
#endif

    default:
      gk_errexit(SIGERR, "Unknown objtype: %d\n", ctrl->objtype);
  }

  gk_free((void **)&ubvec, &bestwhere, LTERM);

}


/*************************************************************************/
/*! Computes the initial k-way partitioning using region growing */
/*************************************************************************/
void InitKWayPartitioningGrow(ctrl_t *ctrl, graph_t *graph)
{
  idx_t ntrials, curobj=0, bestobj=0;
  idx_t *bestwhere=NULL;

  WCOREPUSH;

  bestwhere = iwspacemalloc(ctrl, graph->nvtxs);

  for (ntrials=0; ntrials<ctrl->nIparts; ntrials++) {
    GrowKWayPartitioning(ctrl, graph);

    switch (ctrl->objtype) {
      case METIS_OBJTYPE_CUT:
        curobj = graph->mincut;
        break;

      case METIS_OBJTYPE_VOL:
        curobj = graph->minvol;
        break;

      default:
        gk_errexit(SIGERR, "Unknown objtype: %d\n", ctrl->objtype);
    }
    printf("  Ipart.%"PRIDX" curobj: %"PRIDX"\n", ntrials, curobj);

    if (ntrials == 0 || bestobj > curobj) {
      icopy(graph->nvtxs, graph->where, bestwhere);
      bestobj = curobj;
    }
  }

  if (bestobj != curobj)
    icopy(graph->nvtxs, bestwhere, graph->where);

  printf("  Ipart Select bestobj: %"PRIDX"\n", bestobj);

  WCOREPOP;

}


/*************************************************************************/
/*! Grows and refines a k-way partitioning */
/*************************************************************************/
void GrowKWayPartitioning(ctrl_t *ctrl, graph_t *graph)
{
  idx_t v, i, ii, j, nvtxs, nparts, nmis, minvwgt;
  idx_t *xadj, *adjncy, *vwgt, *adjwgt;
  idx_t *where, *pwgts, *minwgt, *maxwgt, *nbrwgt;
  idx_t *perm;
  real_t ubfactor_original;


  WCOREPUSH;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  vwgt   = graph->vwgt;

  where  = graph->where;
  pwgts  = graph->pwgts;

  nparts = ctrl->nparts;
  
  /* setup the weight intervals of the various subdomains */
  minwgt  = iwspacemalloc(ctrl, nparts);
  maxwgt  = iwspacemalloc(ctrl, nparts);

  for (i=0; i<nparts; i++) {
    maxwgt[i] = ctrl->tpwgts[i]*graph->tvwgt[0]*ctrl->ubfactors[0];
    minwgt[i] = ctrl->tpwgts[i]*graph->tvwgt[0]*(1.0/ctrl->ubfactors[0]);
  }

  /* setup the initial state of the partitioning */
  iset(nvtxs, nparts, where);
  iset(nparts, 0, pwgts);

  perm = iwspacemalloc(ctrl, nvtxs);

  /* compute the weights of the neighborhood */
  nbrwgt = iwspacemalloc(ctrl, nvtxs);
  /*
  for (i=0; i<nvtxs; i++) {
    nbrwgt[i] = vwgt[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      ii = adjncy[j];
      nbrwgt[i] += vwgt[ii]/log2(2+xadj[ii+1]-xadj[ii]);
    }
  }
  */
  for (i=0; i<nvtxs; i++) {
    nbrwgt[i] = 0;
    for (j=xadj[i]; j<xadj[i+1]; j++) 
      nbrwgt[i] += adjwgt[j];
  }
  minvwgt = isum(nvtxs, nbrwgt, 1)/nvtxs;


#ifdef XXX
  perm    = iwspacemalloc(ctrl, nvtxs);
  tperm   = iwspacemalloc(ctrl, nvtxs);
  degrees = iwspacemalloc(ctrl, nvtxs);

  irandArrayPermute(nvtxs, tperm, nvtxs, 1);
  irandArrayPermute(nvtxs, tperm, nvtxs, 0);

  avgdegree = 1.0*(xadj[nvtxs]/nvtxs);
  for (i=0; i<nvtxs; i++) {
    bnum = sqrt(1+xadj[i+1]-xadj[i]);
    degrees[i] = (bnum > avgdegree ? avgdegree : bnum);
  }
  BucketSortKeysInc(ctrl, nvtxs, avgdegree, degrees, tperm, perm);
#endif

  ubfactor_original  = ctrl->ubfactors[0];
  ctrl->ubfactors[0] = 1.05*ubfactor_original;

  /* find an MIS of vertices by randomly traversing the vertices and assign them to
     the different partitions */
  irandArrayPermute(nvtxs, perm, nvtxs, 1);
  for (nmis=0, ii=0; ii<nvtxs; ii++) {
    i=perm[ii];

    if (nbrwgt[i] < minvwgt)
      continue;

    if (where[i] == nparts) {
      pwgts[nmis] = vwgt[i];
      where[i]    = nmis++;

      /* mark first level neighbors */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        v = adjncy[j];
        if (where[v] == nparts)
          where[v] = nparts+1;
      }
    }
    if (nmis == nparts)
      break;
  }
  printf("  nvtxs: %"PRIDX", nmis: %"PRIDX", minvwgt: %"PRIDX", ii: %"PRIDX"\n", nvtxs, nmis, minvwgt, ii);


  /* if the size of the MIS is not sufficiently large, go and find some additional seeds */
  if (nmis < nparts) {
    minvwgt = .75*minvwgt;

    irandArrayPermute(nvtxs, perm, nvtxs, 0);
    for (ii=0; ii<nvtxs; ii++) {
      i=perm[ii];

      if (nbrwgt[i] < minvwgt)
        continue;

      if (where[i] == nparts) {
        pwgts[nmis] = vwgt[i];
        where[i]    = nmis++;

        /* mark first level neighbors */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          v = adjncy[j];
          if (where[v] == nparts)
            where[v] = nparts+1;
        }
      }
      if (nmis == nparts)
        break;
    }

    printf("  nvtxs: %"PRIDX", nmis: %"PRIDX"\n", nvtxs, nmis);
  }

  /* if the size of the MIS is not sufficiently large, go and find some additional seeds */
  if (nmis < nparts) {
    irandArrayPermute(nvtxs, perm, nvtxs, 0);
    for (ii=0; ii<nvtxs; ii++) {
      i = perm[ii];

      if (where[i] == nparts+1) { 
        pwgts[nmis] = vwgt[i];
        where[i]    = nmis++;
      }
      if (nmis == nparts)
        break;
    }
    printf("  nvtxs: %"PRIDX", nmis: %"PRIDX"\n", nvtxs, nmis);
  }

  /* set all unassigned vertices to 'nparts' */
  for (i=0; i<nvtxs; i++) {
    if (where[i] >= nparts)
      where[i] = nparts;
  }

  WCOREPOP;

  /* refine the partition */
  ComputeKWayPartitionParams(ctrl, graph);

  if (ctrl->minconn)
    EliminateSubDomainEdges(ctrl, graph);

  for (i=0; i<4; i++) {
    ComputeKWayBoundary(ctrl, graph, BNDTYPE_BALANCE);
    Greedy_KWayOptimize(ctrl, graph, 10, 0, OMODE_BALANCE);
    /*
    for (k=0; k<nparts; k++)
      printf("%"PRIDX"\n", graph->pwgts[k]);
    exit(0);
    */


    ComputeKWayBoundary(ctrl, graph, BNDTYPE_REFINE);
    Greedy_KWayOptimize(ctrl, graph, ctrl->niter, 1, OMODE_REFINE);
    Greedy_KWayEdgeCutOptimize(ctrl, graph, ctrl->niter);
  }

  ctrl->ubfactors[0] = ubfactor_original;
}
