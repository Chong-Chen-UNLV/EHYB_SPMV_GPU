/*
 * Copyright 1997-2015, Regents of the University of Minnesota
 *
 * rename.h
 *
 * This file contains header files
 *
 * Started 10/2/97
 * George
 *
 *
 * - Dominique LaSalle 2015-02-28
 * Modified to obfuscate names using the 'MTMETIS_' prefix.
 *
 */


#ifndef _LIBMETIS_RENAME_H_
#define _LIBMETIS_RENAME_H_


/* balance.c */
#define Balance2Way			MTMETIS_Balance2Way
#define Bnd2WayBalance			MTMETIS_Bnd2WayBalance
#define General2WayBalance		MTMETIS_General2WayBalance
#define McGeneral2WayBalance            MTMETIS_McGeneral2WayBalance

/* bucketsort.c */
#define BucketSortKeysInc		MTMETIS_BucketSortKeysInc

/* checkgraph.c */
#define CheckGraph                      MTMETIS_CheckGraph
#define CheckInputGraphWeights          MTMETIS_CheckInputGraphWeights
#define FixGraph                        MTMETIS_FixGraph

/* coarsen.c */
#define CoarsenGraph			MTMETIS_CoarsenGraph
#define Match_RM                        MTMETIS_Match_RM
#define Match_SHEM                      MTMETIS_Match_SHEM
#define Match_2Hop                      MTMETIS_Match_2Hop
#define Match_2HopAny                   MTMETIS_Match_2HopAny
#define Match_2HopAll                   MTMETIS_Match_2HopAll
#define PrintCGraphStats                MTMETIS_PrintCGraphStats
#define CreateCoarseGraph		MTMETIS_CreateCoarseGraph
#define CreateCoarseGraphNoMask		MTMETIS_CreateCoarseGraphNoMask
#define CreateCoarseGraphPerm		MTMETIS_CreateCoarseGraphPerm
#define SetupCoarseGraph		MTMETIS_SetupCoarseGraph
#define ReAdjustMemory			MTMETIS_ReAdjustMemory

/* compress.c */
#define CompressGraph			MTMETIS_CompressGraph
#define PruneGraph			MTMETIS_PruneGraph

/* contig.c */
#define FindPartitionInducedComponents  MTMETIS_FindPartitionInducedComponents   
#define IsConnected                     MTMETIS_IsConnected
#define IsConnectedSubdomain            MTMETIS_IsConnectedSubdomain
#define FindSepInducedComponents        MTMETIS_FindSepInducedComponents
#define EliminateComponents             MTMETIS_EliminateComponents
#define MoveGroupContigForCut           MTMETIS_MoveGroupContigForCut
#define MoveGroupContigForVol           MTMETIS_MoveGroupContigForVol

/* debug.c */
#define ComputeCut			MTMETIS_ComputeCut
#define ComputeVolume			MTMETIS_ComputeVolume
#define ComputeMaxCut			MTMETIS_ComputeMaxCut
#define CheckBnd			MTMETIS_CheckBnd
#define CheckBnd2			MTMETIS_CheckBnd2
#define CheckNodeBnd			MTMETIS_CheckNodeBnd
#define CheckRInfo			MTMETIS_CheckRInfo
#define CheckNodePartitionParams	MTMETIS_CheckNodePartitionParams
#define IsSeparable			MTMETIS_IsSeparable
#define CheckKWayVolPartitionParams     MTMETIS_CheckKWayVolPartitionParams

/* fm.c */
#define FM_2WayRefine                   MTMETIS_FM_2WayRefine
#define FM_2WayCutRefine                MTMETIS_FM_2WayCutRefine
#define FM_Mc2WayCutRefine              MTMETIS_FM_Mc2WayCutRefine
#define SelectQueue                     MTMETIS_SelectQueue
#define Print2WayRefineStats            MTMETIS_Print2WayRefineStats

/* fortran.c */
#define Change2CNumbering		MTMETIS_Change2CNumbering
#define Change2FNumbering		MTMETIS_Change2FNumbering
#define Change2FNumbering2		MTMETIS_Change2FNumbering2
#define Change2FNumberingOrder		MTMETIS_Change2FNumberingOrder
#define ChangeMesh2CNumbering		MTMETIS_ChangeMesh2CNumbering
#define ChangeMesh2FNumbering		MTMETIS_ChangeMesh2FNumbering
#define ChangeMesh2FNumbering2		MTMETIS_ChangeMesh2FNumbering2

/* graph.c */
#define SetupGraph			MTMETIS_SetupGraph
#define SetupGraph_adjrsum              MTMETIS_SetupGraph_adjrsum
#define SetupGraph_tvwgt                MTMETIS_SetupGraph_tvwgt
#define SetupGraph_label                MTMETIS_SetupGraph_label
#define SetupSplitGraph                 MTMETIS_SetupSplitGraph
#define CreateGraph                     MTMETIS_CreateGraph
#define InitGraph                       MTMETIS_InitGraph
#define FreeRData                       MTMETIS_FreeRData
#define FreeGraph                       MTMETIS_FreeGraph
#define graph_WriteToDisk               MTMETIS_graph_WriteToDisk
#define graph_ReadFromDisk              MTMETIS_graph_ReadFromDisk

/* initpart.c */
#define Init2WayPartition		MTMETIS_Init2WayPartition
#define InitSeparator			MTMETIS_InitSeparator
#define RandomBisection			MTMETIS_RandomBisection
#define GrowBisection			MTMETIS_GrowBisection
#define McRandomBisection               MTMETIS_McRandomBisection
#define McGrowBisection                 MTMETIS_McGrowBisection
#define GrowBisectionNode		MTMETIS_GrowBisectionNode

/* kmetis.c */
#define MlevelKWayPartitioning		MTMETIS_MlevelKWayPartitioning
#define InitKWayPartitioning            MTMETIS_InitKWayPartitioning
#define InitKWayPartitioningRB          MTMETIS_InitKWayPartitioningRB
#define InitKWayPartitioningGrow        MTMETIS_InitKWayPartitioningGrow

/* kwayfm.c */
#define Greedy_KWayOptimize		MTMETIS_Greedy_KWayOptimize
#define Greedy_KWayCutOptimize		MTMETIS_Greedy_KWayCutOptimize
#define Greedy_KWayVolOptimize          MTMETIS_Greedy_KWayVolOptimize
#define Greedy_McKWayCutOptimize        MTMETIS_Greedy_McKWayCutOptimize
#define Greedy_McKWayVolOptimize        MTMETIS_Greedy_McKWayVolOptimize
#define IsArticulationNode              MTMETIS_IsArticulationNode
#define KWayVolUpdate                   MTMETIS_KWayVolUpdate

/* kwayrefine.c */
#define RefineKWay			MTMETIS_RefineKWay
#define AllocateKWayPartitionMemory	MTMETIS_AllocateKWayPartitionMemory
#define ComputeKWayPartitionParams	MTMETIS_ComputeKWayPartitionParams
#define ProjectKWayPartition		MTMETIS_ProjectKWayPartition
#define ComputeKWayBoundary		MTMETIS_ComputeKWayBoundary
#define ComputeKWayVolGains             MTMETIS_ComputeKWayVolGains
#define IsBalanced			MTMETIS_IsBalanced

/* mcutil */
#define rvecle                          MTMETIS_rvecle
#define rvecge                          MTMETIS_rvecge
#define rvecsumle                       MTMETIS_rvecsumle
#define rvecmaxdiff                     MTMETIS_rvecmaxdiff
#define ivecle                          MTMETIS_ivecle
#define ivecge                          MTMETIS_ivecge
#define ivecaxpylez                     MTMETIS_ivecaxpylez
#define ivecaxpygez                     MTMETIS_ivecaxpygez
#define BetterVBalance                  MTMETIS_BetterVBalance
#define BetterBalance2Way               MTMETIS_BetterBalance2Way
#define BetterBalanceKWay               MTMETIS_BetterBalanceKWay
#define ComputeLoadImbalance            MTMETIS_ComputeLoadImbalance
#define ComputeLoadImbalanceDiff        MTMETIS_ComputeLoadImbalanceDiff
#define ComputeLoadImbalanceDiffVec     MTMETIS_ComputeLoadImbalanceDiffVec
#define ComputeLoadImbalanceVec         MTMETIS_ComputeLoadImbalanceVec

/* mesh.c */
#define CreateGraphDual                 MTMETIS_CreateGraphDual
#define FindCommonElements              MTMETIS_FindCommonElements
#define CreateGraphNodal                MTMETIS_CreateGraphNodal
#define FindCommonNodes                 MTMETIS_FindCommonNodes
#define CreateMesh                      MTMETIS_CreateMesh
#define InitMesh                        MTMETIS_InitMesh
#define FreeMesh                        MTMETIS_FreeMesh

/* meshpart.c */
#define InduceRowPartFromColumnPart     MTMETIS_InduceRowPartFromColumnPart

/* minconn.c */
#define ComputeSubDomainGraph           MTMETIS_ComputeSubDomainGraph
#define UpdateEdgeSubDomainGraph        MTMETIS_UpdateEdgeSubDomainGraph
#define PrintSubDomainGraph             MTMETIS_PrintSubDomainGraph
#define EliminateSubDomainEdges         MTMETIS_EliminateSubDomainEdges
#define MoveGroupMinConnForCut          MTMETIS_MoveGroupMinConnForCut
#define MoveGroupMinConnForVol          MTMETIS_MoveGroupMinConnForVol

/* mincover.c */
#define MinCover			MTMETIS_MinCover
#define MinCover_Augment		MTMETIS_MinCover_Augment
#define MinCover_Decompose		MTMETIS_MinCover_Decompose
#define MinCover_ColDFS			MTMETIS_MinCover_ColDFS
#define MinCover_RowDFS			MTMETIS_MinCover_RowDFS

/* mmd.c */
#define genmmd				MTMETIS_genmmd
#define mmdelm				MTMETIS_mmdelm
#define mmdint				MTMETIS_mmdint
#define mmdnum				MTMETIS_mmdnum
#define mmdupd				MTMETIS_mmdupd


/* ometis.c */
#define MlevelNestedDissection		MTMETIS_MlevelNestedDissection
#define MlevelNestedDissectionCC	MTMETIS_MlevelNestedDissectionCC
#define MlevelNodeBisectionMultiple	MTMETIS_MlevelNodeBisectionMultiple
#define MlevelNodeBisectionL2		MTMETIS_MlevelNodeBisectionL2
#define MlevelNodeBisectionL1		MTMETIS_MlevelNodeBisectionL1
#define SplitGraphOrder			MTMETIS_SplitGraphOrder
#define SplitGraphOrderCC		MTMETIS_SplitGraphOrderCC
#define MMDOrder			MTMETIS_MMDOrder

/* options.c */
#define SetupCtrl                       MTMETIS_SetupCtrl
#define SetupKWayBalMultipliers         MTMETIS_SetupKWayBalMultipliers
#define Setup2WayBalMultipliers         MTMETIS_Setup2WayBalMultipliers
#define PrintCtrl                       MTMETIS_PrintCtrl
#define FreeCtrl                        MTMETIS_FreeCtrl
#define CheckParams                     MTMETIS_CheckParams

/* parmetis.c */
#define MlevelNestedDissectionP		MTMETIS_MlevelNestedDissectionP
#define FM_2WayNodeRefine1SidedP        MTMETIS_FM_2WayNodeRefine1SidedP
#define FM_2WayNodeRefine2SidedP        MTMETIS_FM_2WayNodeRefine2SidedP

/* pmetis.c */
#define MlevelRecursiveBisection	MTMETIS_MlevelRecursiveBisection
#define MultilevelBisect		MTMETIS_MultilevelBisect
#define SplitGraphPart			MTMETIS_SplitGraphPart

/* refine.c */
#define Refine2Way			MTMETIS_Refine2Way
#define Allocate2WayPartitionMemory	MTMETIS_Allocate2WayPartitionMemory
#define Compute2WayPartitionParams	MTMETIS_Compute2WayPartitionParams
#define Project2WayPartition		MTMETIS_Project2WayPartition

/* separator.c */
#define ConstructSeparator		MTMETIS_ConstructSeparator
#define ConstructMinCoverSeparator	MTMETIS_ConstructMinCoverSeparator

/* sfm.c */
#define FM_2WayNodeRefine2Sided         MTMETIS_FM_2WayNodeRefine2Sided 
#define FM_2WayNodeRefine1Sided         MTMETIS_FM_2WayNodeRefine1Sided
#define FM_2WayNodeBalance              MTMETIS_FM_2WayNodeBalance

/* srefine.c */
#define Refine2WayNode			MTMETIS_Refine2WayNode
#define Allocate2WayNodePartitionMemory	MTMETIS_Allocate2WayNodePartitionMemory
#define Compute2WayNodePartitionParams	MTMETIS_Compute2WayNodePartitionParams
#define Project2WayNodePartition	MTMETIS_Project2WayNodePartition

/* stat.c */
#define ComputePartitionInfoBipartite   MTMETIS_ComputePartitionInfoBipartite
#define ComputePartitionBalance		MTMETIS_ComputePartitionBalance
#define ComputeElementBalance		MTMETIS_ComputeElementBalance

/* timing.c */
#define InitTimers			MTMETIS_InitTimers
#define PrintTimers			MTMETIS_PrintTimers

/* util.c */
#define irandInRange                    MTMETIS_irandInRange
#define irandArrayPermute               MTMETIS_irandArrayPermute
#define iargmax_strd                    MTMETIS_iargmax_strd 
#define iargmax_nrm                     MTMETIS_iargmax_nrm
#define iargmax2_nrm                    MTMETIS_iargmax2_nrm
#define rargmax2                        MTMETIS_rargmax2
#define InitRandom                      MTMETIS_InitRandom
#define metis_rcode                     MTMETIS_metis_rcode

/* wspace.c */
#define AllocateWorkSpace               MTMETIS_AllocateWorkSpace                  
#define AllocateRefinementWorkSpace     MTMETIS_AllocateRefinementWorkSpace
#define FreeWorkSpace                   MTMETIS_FreeWorkSpace
#define wspacemalloc                    MTMETIS_wspacemalloc
#define wspacepush                      MTMETIS_wspacepush
#define wspacepop                       MTMETIS_wspacepop
#define iwspacemalloc                   MTMETIS_iwspacemalloc
#define rwspacemalloc                   MTMETIS_rwspacemalloc
#define ikvwspacemalloc                 MTMETIS_ikvwspacemalloc
#define cnbrpoolReset                   MTMETIS_cnbrpoolReset
#define cnbrpoolGetNext                 MTMETIS_cnbrpoolGetNext
#define vnbrpoolReset                   MTMETIS_vnbrpoolReset
#define vnbrpoolGetNext                 MTMETIS_vnbrpoolGetNext


/* Hide toplevel names */
#define METIS_PartGraphRecursive __METIS_PartGraphRecursive 
#define METIS_PartGraphKway __METIS_PartGraphKway
#define METIS_NodeND __METIS_NodeND
#define METIS_Free __METIS_Free
#define METIS_SetDefaultOptions __METIS_SetDefaultOptions
#define METIS_ComputeVertexSeparator __METIS_ComputeVertexSeparator
#define METIS_NodeRefine __METIS_NodeRefine
 

#endif


