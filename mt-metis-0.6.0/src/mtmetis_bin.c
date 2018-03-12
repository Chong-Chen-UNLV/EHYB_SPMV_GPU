/**
 * @file mtmetis_bin.c
 * @brief Main driver function
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2014, Regents of the University of Minnesota
 * @version 1
 * @date 2013-05-20
 */





#ifndef MTMETIS_BIN_C
#define MTMETIS_BIN_C




#include "base.h"
#include "partition.h"
#include "order.h"
#include "graph.h"
#include "ctrl.h"

#include <wildriver.h>



/******************************************************************************
* MACROS **********************************************************************
******************************************************************************/


#define S_ARRAY_SIZE(a) \
  (sizeof(a) > 0 ? (sizeof(a) / sizeof((a)[0])) : 0)




/******************************************************************************
* OPTIONS *********************************************************************
******************************************************************************/


static const cmd_opt_pair_t CTYPE_CHOICES[] = {
  {MTMETIS_STR_CTYPE_RM,"Random Matching",MTMETIS_CTYPE_RM},
  {MTMETIS_STR_CTYPE_SHEM,"Sorted Heavy Edge Matching",MTMETIS_CTYPE_SHEM},
  {MTMETIS_STR_CTYPE_FC,"FirstChoice Grouping",MTMETIS_CTYPE_FC}
};


static const cmd_opt_pair_t CONTYPE_CHOICES[] = {
  {MTMETIS_STR_CONTYPE_CLS,"Hash-table with linear scanning", \
      MTMETIS_CONTYPE_CLS},
  {MTMETIS_STR_CONTYPE_DENSE,"Dense vector",MTMETIS_CONTYPE_DENSE},
  {MTMETIS_STR_CONTYPE_SORT,"Sort and merge",MTMETIS_CONTYPE_SORT}
};


static const cmd_opt_pair_t RTYPE_CHOICES[] = {
  {MTMETIS_STR_RTYPE_GREEDY,"Greedy refinement",MTMETIS_RTYPE_GREEDY},
  {MTMETIS_STR_RTYPE_FM,"FM based serial refinement (pairwise for kway)", \
      MTMETIS_RTYPE_FM},
  {MTMETIS_STR_RTYPE_SFM,"Segmented FM based parallel refinement", \
      MTMETIS_RTYPE_SFM},
  {MTMETIS_STR_RTYPE_SFG,"Segmented FM plus Greedy parallel refinement", \
      MTMETIS_RTYPE_SFG},
  {MTMETIS_STR_RTYPE_HS,"Hill-Scanning refinement", \
      MTMETIS_RTYPE_HS}
};


static const cmd_opt_pair_t PTYPE_CHOICES[] = {
  {MTMETIS_STR_PTYPE_KWAY,"K-Way Edgecut",MTMETIS_PTYPE_KWAY},
  {MTMETIS_STR_PTYPE_ESEP,"Edge Separator",MTMETIS_PTYPE_ESEP},
  {MTMETIS_STR_PTYPE_RB,"RB Edgecut",MTMETIS_PTYPE_RB},
  {MTMETIS_STR_PTYPE_VSEP,"Vertex Separator",MTMETIS_PTYPE_VSEP},
  {MTMETIS_STR_PTYPE_ND,"Nested Dissection",MTMETIS_PTYPE_ND}
};


static const cmd_opt_pair_t VERBOSITY_CHOICES[] = {
  {MTMETIS_STR_VERBOSITY_NONE,"Do not print any runtime information.", \
      MTMETIS_VERBOSITY_NONE},
  {MTMETIS_STR_VERBOSITY_LOW,"Print only metric summary.", \
      MTMETIS_VERBOSITY_LOW},
  {MTMETIS_STR_VERBOSITY_MEDIUM,"Print summary run information.", \
      MTMETIS_VERBOSITY_MEDIUM},
  {MTMETIS_STR_VERBOSITY_HIGH,"Print verbose run information.", \
      MTMETIS_VERBOSITY_HIGH},
  {MTMETIS_STR_VERBOSITY_MAXIMUM,"Print everything.", \
      MTMETIS_VERBOSITY_MAXIMUM}
};


static const cmd_opt_pair_t DISTRIBUTION_CHOICES[] = {
  {MTMETIS_STR_DISTRIBUTION_BLOCK,"Distribute the vertices in continuous " \
      "from the initial ordering.",MTMETIS_DISTRIBUTION_BLOCK},
  {MTMETIS_STR_DISTRIBUTION_CYCLIC,"Distribute the vertices in a cyclic " \
      "fashion.",MTMETIS_DISTRIBUTION_CYCLIC},
  {MTMETIS_STR_DISTRIBUTION_BLOCKCYCLIC,"Distribute the vertices in a " \
      "blockcyclic fashion.",MTMETIS_DISTRIBUTION_BLOCKCYCLIC}
};


static const cmd_opt_pair_t IGNOREWEIGHTS_CHOICES[] = {
  {MTMETIS_STR_IGNORE_NONE,"Use all weights normally", \
      MTMETIS_IGNORE_NONE},
  {MTMETIS_STR_IGNORE_VERTEXWEIGHTS,"Force all vertex weights to be one", \
      MTMETIS_IGNORE_VERTEXWEIGHTS},
  {MTMETIS_STR_IGNORE_EDGEWEIGHTS,"Force all edge weights to be one", \
      MTMETIS_IGNORE_EDGEWEIGHTS},
  {MTMETIS_STR_IGNORE_BOTH,"Force all weights to be one", \
      MTMETIS_IGNORE_VERTEXWEIGHTS | MTMETIS_IGNORE_EDGEWEIGHTS}
};


static const cmd_opt_pair_t SCANTYPE_CHOICES[] = {
  {MTMETIS_STR_SCANTYPE_SQRT,"Use the square-root of the number of boundary " \
      "vertices.",MTMETIS_HS_SCAN_SQRT},
  {MTMETIS_STR_SCANTYPE_1PC,"Use 1% of the number of boundary " \
      "vertices.",MTMETIS_HS_SCAN_1PC},
  {MTMETIS_STR_SCANTYPE_5PC,"Use 5% of the number of boundary " \
      "vertices.",MTMETIS_HS_SCAN_5PC},
  {MTMETIS_STR_SCANTYPE_25PC,"Use 25% of the number of boundary " \
      "vertices.",MTMETIS_HS_SCAN_25PC},
  {MTMETIS_STR_SCANTYPE_FULL,"Do a full-scan, no limit.",MTMETIS_HS_SCAN_FULL},
};


static const cmd_opt_t OPTS[] = {
  {MTMETIS_OPTION_HELP,'h',"help","Display this help page.",CMD_OPT_FLAG, \
      NULL,0},
  {MTMETIS_OPTION_CTYPE,'c',"ctype","The type of coarsening (default=shem).", \
      CMD_OPT_CHOICE,CTYPE_CHOICES,S_ARRAY_SIZE(CTYPE_CHOICES)},
  {MTMETIS_OPTION_CONTYPE,'d',"contype","How to merge adjacency lists " \
      "during contraction (default=ls).",CMD_OPT_CHOICE,CONTYPE_CHOICES, \
        S_ARRAY_SIZE(CONTYPE_CHOICES)},
  {MTMETIS_OPTION_RTYPE,'r',"rtype","The type of refinement " \
      "(default=greedy).",CMD_OPT_CHOICE,RTYPE_CHOICES, \
      S_ARRAY_SIZE(RTYPE_CHOICES)},
  {MTMETIS_OPTION_SEED,'s',"seed","The random seed to use.",CMD_OPT_INT,NULL, \
      0},
  {MTMETIS_OPTION_NCUTS,'N',"cuts","The number of cuts to " \
      "generate using successive random seeds (default=1).",CMD_OPT_INT,NULL, \
      0},
  {MTMETIS_OPTION_NRUNS,'n',"runs","The number of partitionings to " \
      "generate using successive random seeds (default=1).",CMD_OPT_INT,NULL, \
      0},
  {MTMETIS_OPTION_NINITSOLUTIONS,'i',"initialcuts","The number of " \
      "initial cuts to generate at the coarsest level (default=8).", \
      CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_NITER,'R',"nrefpass","The maximum number of refinement " \
      "passes (default=8).",CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_TIME,'t',"times","Print timing information",CMD_OPT_FLAG, \
      NULL,0},
  {MTMETIS_OPTION_NTHREADS,'T',"threads","The number of threads to use.", \
      CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_VERBOSITY,'v',"verbosity","The amount of information to " \
      "print during partitioning (default=none).",CMD_OPT_CHOICE, \
      VERBOSITY_CHOICES,S_ARRAY_SIZE(VERBOSITY_CHOICES)},
  {MTMETIS_OPTION_DISTRIBUTION,'D',"distribution","The distribution to use " \
      "for assigning vertices to threads (default=blockcyclic).", \
      CMD_OPT_CHOICE,DISTRIBUTION_CHOICES,S_ARRAY_SIZE(DISTRIBUTION_CHOICES)},
  {MTMETIS_OPTION_UBFACTOR,'b',"balance","The balance constraint " \
      "(default=1.03, which means allowing for a 3% imbalance).", \
      CMD_OPT_FLOAT,NULL,0},
  {MTMETIS_OPTION_PTYPE,'p',"ptype","The type of partition to compute " \
      "(default=kway)",CMD_OPT_CHOICE,PTYPE_CHOICES, \
      S_ARRAY_SIZE(PTYPE_CHOICES)},
  {MTMETIS_OPTION_RUNSTATS,'C',"partstats","Print statics on quality of " \
      "partitions.",CMD_OPT_FLAG,NULL,0},
  {MTMETIS_OPTION_METIS,'M',"metis","When run with one thread, call Metis " \
      "directly.",CMD_OPT_FLAG,NULL,0},
  {MTMETIS_OPTION_LEAFMATCH,'L',"leafmatch","Match leaf vertices together " \
      "if there are too many unmatched vertices (default=true).", \
      CMD_OPT_BOOL,NULL,0},
  {MTMETIS_OPTION_REMOVEISLANDS,'I',"removeislands","Remove island vertices " \
      "before partitioning (default=false).",CMD_OPT_BOOL,NULL,0},
  {MTMETIS_OPTION_VWGTDEGREE,'V',"vwgtdegree","Use the degree of each " \
      "vertex as its weight (default=false).",CMD_OPT_FLAG,NULL,0},
  {MTMETIS_OPTION_IGNORE,'W',"ignoreweights","Ignore input weights " \
      "on a graph file (default=none).",CMD_OPT_CHOICE,IGNOREWEIGHTS_CHOICES, \
      S_ARRAY_SIZE(IGNOREWEIGHTS_CHOICES)},
  {MTMETIS_OPTION_HILLSIZE,'H',"hillsize","The limit to use when searching " \
      "for hills (default=16). This only applies to hill climbing " \
      "refinement types.",CMD_OPT_INT,NULL,0},
  {MTMETIS_OPTION_HS_SCANTYPE,'S',"scantype","The how many hills to scan " \
      "before terminating (default=sqrt). This only applies to HS " \
      "refinement.",CMD_OPT_CHOICE, \
      SCANTYPE_CHOICES,S_ARRAY_SIZE(SCANTYPE_CHOICES)},
  {MTMETIS_OPTION_VERSION,'\0',"version","Display the current version.", \
      CMD_OPT_FLAG,NULL,0}
};


static const size_t NOPTS = S_ARRAY_SIZE(OPTS);


#undef S_ARRAY_SIZE




/******************************************************************************
* PRIVATE FUNCTIONS ***********************************************************
******************************************************************************/


static int S_usage(
    char const * const name,
    FILE * const fout)
{
  fprintf(fout,"USAGE:\n");
  fprintf(fout,"%s [options] <graphfile> <nparts> [ <partfile> | - ]\n", \
      name);
  fprintf(fout,"%s -p nd [options] <graphfile> [ <permfile> | - ]\n", \
      name);
  fprintf(fout,"\n");
  fprintf(fout,"Options:\n");
  fprint_cmd_opts(fout,OPTS,NOPTS);

  return 1;
}


static double * S_parse_args(
    cmd_arg_t * args, 
    size_t nargs,
    char const ** r_input, 
    char const ** r_output)
{
  size_t i, xarg;
  double * options = NULL;
  const char * input_file = NULL, * output_file = NULL;

  options = mtmetis_init_options();

  /* set default verbosity to low */
  options[MTMETIS_OPTION_VERBOSITY] = MTMETIS_VERBOSITY_LOW;
  options[MTMETIS_OPTION_NPARTS] = 2.0;
  options[MTMETIS_OPTION_PTYPE] = MTMETIS_PTYPE_KWAY;

  for (i=0;i<nargs;++i) {
    switch (args[i].type) {
      case CMD_OPT_CHOICE:
        options[args[i].id] = (double)args[i].val.o;
        break;
      case CMD_OPT_BOOL:
        options[args[i].id] = (double)args[i].val.b;
        break;
      case CMD_OPT_INT:
        options[args[i].id] = (double)args[i].val.i;
        break;
      case CMD_OPT_FLOAT:
        options[args[i].id] = (double)args[i].val.f;
        break;
      case CMD_OPT_FLAG:
        options[args[i].id] = 1.0;
        break;
      default:
        break;
    }
  }

  xarg = 0;
  for (i=0;i<nargs;++i) {
    /* check for help */
    if (args[i].id == MTMETIS_OPTION_HELP) {
      goto CLEANUP;
    } else if (args[i].id == MTMETIS_OPTION_VERSION) {
      printf("mt-Metis %d.%d.%d\n",MTMETIS_VER_MAJOR,MTMETIS_VER_MINOR, \
          MTMETIS_VER_SUBMINOR);
      printf("Copyright 2016, The Regents of the University of Minnesota\n");
      goto CLEANUP;
    }
    if (args[i].type == CMD_OPT_XARG) {
      if (xarg == 0) {
        input_file = args[i].val.s;
      } else {
        if (options[MTMETIS_OPTION_PTYPE] == MTMETIS_PTYPE_KWAY || \
            options[MTMETIS_OPTION_PTYPE] == MTMETIS_PTYPE_RB) {
          if (xarg == 1) {
            options[MTMETIS_OPTION_NPARTS] = (pid_type)atoll(args[i].val.s);
          } else if (xarg == 2) {
            output_file = args[i].val.s;
            if (strcmp(output_file,"-") == 0) {
              /* if we are going to print to stdout, don't print anything else */
              options[MTMETIS_OPTION_VERBOSITY] = MTMETIS_VERBOSITY_NONE;
            }
          } else {
            eprintf("Unknown extra argument '%s'\n",args[i].val.s);
            goto CLEANUP;
          }
        } else {
          if (xarg == 1) {
            output_file = args[i].val.s;
            if (strcmp(output_file,"-") == 0) {
              /* if we are going to print to stdout, don't print anything else */
              options[MTMETIS_OPTION_VERBOSITY] = MTMETIS_VERBOSITY_NONE;
            }
          } else {
            eprintf("Unknown extra argument '%s'\n",args[i].val.s);
            goto CLEANUP;
          }
        }
      }
      ++xarg;
    }
  }

  if (input_file == NULL) {
    eprintf("Must supply at least an input graph to partition\n");
    goto CLEANUP;
  }

  *r_output = output_file;
  *r_input = input_file;

  return options;

  CLEANUP:
  dl_free(options);
  *r_output = NULL;
  *r_input = NULL;

  return NULL;
}




/******************************************************************************
* MAIN ************************************************************************
******************************************************************************/


int main(
    int argc, 
    char ** argv) 
{
  int rv, times, verbosity;
  size_t nargs;
  vtx_type nvtxs, i;
  adj_type * xadj = NULL;
  vtx_type * adjncy = NULL;
  wgt_type * vwgt = NULL, * adjwgt = NULL;
  double * options = NULL;
  cmd_arg_t * args = NULL;
  pid_type * owhere = NULL;
  char const * output_file = NULL, * input_file = NULL;
  dl_timer_t timer_input, timer_output;

  /* parse user specified options */
  rv = cmd_parse_args(argc-1,argv+1,OPTS,NOPTS,&args,&nargs);
  if (rv != DL_CMDLINE_SUCCESS) {
    S_usage(argv[0],stderr);
    rv = 1;
    goto CLEANUP;
  }
  options = S_parse_args(args,nargs,&input_file,&output_file);
  if (options == NULL) {
    S_usage(argv[0],stderr);
    rv = 2;
    goto CLEANUP;
  }

  /* parse verbosity and timing */
  times = options[MTMETIS_OPTION_TIME];
  verbosity = options[MTMETIS_OPTION_VERBOSITY];

  dl_init_timer(&timer_input);
  dl_init_timer(&timer_output);

  /* start timers */
  dl_start_timer(&timer_input);

  /* read the input graph */
  rv = wildriver_read_graph(input_file,&nvtxs,NULL,NULL,NULL,&xadj,&adjncy, \
      &vwgt,&adjwgt);

  if (rv != 1) {
    rv = 4;
    goto CLEANUP;
  }

  vprintf(verbosity,MTMETIS_VERBOSITY_LOW,"Read '%s' with %"PF_VTX_T \
      " vertices and %"PF_ADJ_T" edges.\n",input_file,nvtxs,xadj[nvtxs]/2);

  dl_stop_timer(&timer_input);

  if (output_file) {
    owhere = pid_alloc(nvtxs);
  }

  if (mtmetis_partition_explicit(nvtxs,xadj,adjncy,vwgt,adjwgt,options,
      owhere,NULL) != MTMETIS_SUCCESS) {
    rv = 3;
    goto CLEANUP;
  }

  dl_start_timer(&timer_output);

  if (output_file) {
    if (strcmp(output_file,"-") == 0) {
      /* write to stdout */
      for (i=0;i<nvtxs;++i) {
        printf("%"PF_PID_T"\n",owhere[i]);
      }
    } else {
      /* save to file */
      FILE * fout = fopen(output_file,"w");
      for (i=0;i<nvtxs;++i) {
        fprintf(fout,"%"PF_PID_T"\n",owhere[i]);
      }
      fclose(fout);
    }
  }

  dl_stop_timer(&timer_output);

  if (times) {
    dl_print_header("AUXILLARY TIME",'$');
    printf("Input: %.05fs\n",dl_poll_timer(&timer_input));
    printf("Output: %.05fs\n",dl_poll_timer(&timer_output));
    dl_print_footer('$');
  }

  CLEANUP:

  if (options) {
    dl_free(options);
  }
  if (xadj) {
    dl_free(xadj);
  }
  if (adjncy) {
    dl_free(adjncy);
  }
  if (vwgt) {
    dl_free(vwgt);
  }
  if (adjwgt) {
    dl_free(adjwgt);
  }
  if (owhere) {
    dl_free(owhere);
  }
  if (args) {
    dl_free(args);
  }

  return 0;
}




#endif
