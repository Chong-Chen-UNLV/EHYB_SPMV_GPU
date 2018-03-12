/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * timing.c
 *
 * This file contains routines that deal with timing Metis
 *
 * Started 7/24/97
 * George
 *
 * $Id: timing.c 12542 2012-08-23 04:49:01Z dominique $
 *
 */

#include "metislib.h"


/*************************************************************************
* This function clears the timers
**************************************************************************/
void InitTimers(ctrl_t *ctrl)
{
  gk_clearwctimer(ctrl->TotalTmr);
  gk_clearwctimer(ctrl->InitPartTmr);
  gk_clearwctimer(ctrl->MatchTmr);
  gk_clearwctimer(ctrl->ContractTmr);
  gk_clearwctimer(ctrl->CoarsenTmr);
  gk_clearwctimer(ctrl->UncoarsenTmr);
  gk_clearwctimer(ctrl->RefTmr);
  gk_clearwctimer(ctrl->ProjectTmr);
  gk_clearwctimer(ctrl->SplitTmr);
  gk_clearwctimer(ctrl->Aux1Tmr);
  gk_clearwctimer(ctrl->Aux2Tmr);
  gk_clearwctimer(ctrl->Aux3Tmr);
}



/*************************************************************************
* This function prints the various timers
**************************************************************************/
void PrintTimers(ctrl_t *ctrl)
{
  printf("\nTiming Information -------------------------------------------------");
  printf("\n Multilevel: \t\t %7.3"PRREAL"", gk_getwctimer(ctrl->TotalTmr));
  printf("\n     Coarsening: \t\t %7.3"PRREAL"", gk_getwctimer(ctrl->CoarsenTmr));
  printf("\n            Matching: \t\t\t %7.3"PRREAL"", gk_getwctimer(ctrl->MatchTmr));
  printf("\n            Contract: \t\t\t %7.3"PRREAL"", gk_getwctimer(ctrl->ContractTmr));
  printf("\n     Initial Partition: \t %7.3"PRREAL"", gk_getwctimer(ctrl->InitPartTmr));
  printf("\n     Uncoarsening: \t\t %7.3"PRREAL"", gk_getwctimer(ctrl->UncoarsenTmr));
  printf("\n          Refinement: \t\t\t %7.3"PRREAL"", gk_getwctimer(ctrl->RefTmr));
  printf("\n          Projection: \t\t\t %7.3"PRREAL"", gk_getwctimer(ctrl->ProjectTmr));
  printf("\n     Splitting: \t\t %7.3"PRREAL"", gk_getwctimer(ctrl->SplitTmr));
/*
  printf("\n       Aux1Tmr: \t\t %7.3"PRREAL"", gk_getwctimer(ctrl->Aux1Tmr));
  printf("\n       Aux2Tmr: \t\t %7.3"PRREAL"", gk_getwctimer(ctrl->Aux2Tmr));
  printf("\n       Aux3Tmr: \t\t %7.3"PRREAL"", gk_getwctimer(ctrl->Aux3Tmr));
*/
  printf("\n********************************************************************\n");
}



