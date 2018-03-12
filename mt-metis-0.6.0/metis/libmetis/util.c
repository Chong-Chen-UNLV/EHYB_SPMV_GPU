/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * util.c
 *
 * This function contains various utility routines
 *
 * Started 9/28/95
 * George
 *
 * $Id: util.c 17622 2014-09-09 03:27:49Z dominique $
 *
 * Dominique LaSalle 2015-02-28
 * Modified to use rentrant random functions.
 *
 */

#define _XOPEN_SOURCE
#include <stdlib.h>

#include "metislib.h"

idx_t my_irandInRange_r(size_t range, idx_t * seed)
{
  uint16_t state[3];
  state[0] = (uint16_t)(*seed);
  state[1] = (uint16_t)((*seed) >> 16);
  state[2] = (uint16_t)((*seed) >> 8);

  uint32_t const r = ((uint32_t)nrand48(state));

  *seed = state[0] + (state[1] << 16);

  return r % range;
}

idx_t my_irandArrayPermute_r(size_t n, idx_t * p, size_t nshuffles,
    idx_t flag, idx_t * seed)
{
  size_t i, u, v;
  idx_t tmp;

  if (flag == 1) {
    for (i=0; i<n; i++)
      p[i] = i;
  }

  if (n < 10) {
    for (i=0; i<n; i++) {
      v = my_irandInRange_r(n,seed);
      u = my_irandInRange_r(n,seed);
      gk_SWAP(p[v], p[u], tmp);
    }
  } else {
    for (i=0; i<nshuffles; i++) {
      v = my_irandInRange_r(n-3,seed);
      u = my_irandInRange_r(n-3,seed);
      gk_SWAP(p[v+0], p[u+2], tmp);
      gk_SWAP(p[v+1], p[u+3], tmp);
      gk_SWAP(p[v+2], p[u+0], tmp);
      gk_SWAP(p[v+3], p[u+1], tmp);
    }
  }

  return 1;
}

/*************************************************************************/
/*! This function initializes the random number generator 
  */
/*************************************************************************/
void InitRandom(idx_t seed)
{
  isrand((seed == -1 ? 4321 : seed)); 
}


/*************************************************************************/
/*! Returns the highest weight index of x[i]*y[i] 
 */
/*************************************************************************/
idx_t iargmax_nrm(size_t n, idx_t *x, real_t *y)
{
  idx_t i, max=0;
      
  for (i=1; i<(idx_t)n; i++)
     max = (x[i]*y[i] > x[max]*y[max] ? i : max);
                
  return max;
}


/*************************************************************************/
/*! These functions return the index of the maximum element in a vector
  */
/*************************************************************************/
idx_t iargmax_strd(size_t n, idx_t *x, idx_t incx)
{
  size_t i, max=0;

  n *= incx;
  for (i=incx; i<n; i+=incx)
    max = (x[i] > x[max] ? i : max);

  return max/incx;
}


/*************************************************************************/
/*! These functions return the index of the almost maximum element in a 
    vector
 */
/*************************************************************************/
idx_t rargmax2(size_t n, real_t *x)
{
  size_t i, max1, max2;

  if (x[0] > x[1]) {
    max1 = 0;
    max2 = 1;
  }
  else {
    max1 = 1;
    max2 = 0;
  }

  for (i=2; i<n; i++) {
    if (x[i] > x[max1]) {
      max2 = max1;
      max1 = i;
    }
    else if (x[i] > x[max2])
      max2 = i;
  }

  return max2;
}


/*************************************************************************/
/*! These functions return the index of the second largest elements in the
    vector formed by x.y where '.' is element-wise multiplication */
/*************************************************************************/
idx_t iargmax2_nrm(size_t n, idx_t *x, real_t *y)
{
  size_t i, max1, max2;

  if (x[0]*y[0] > x[1]*y[1]) {
    max1 = 0;
    max2 = 1;
  }
  else {
    max1 = 1;
    max2 = 0;
  }

  for (i=2; i<n; i++) {
    if (x[i]*y[i] > x[max1]*y[max1]) {
      max2 = max1;
      max1 = i;
    }
    else if (x[i]*y[i] > x[max2]*y[max2])
      max2 = i;
  }

  return max2;
}


/*************************************************************************/
/*! converts a signal code into a Metis return code 
 */
/*************************************************************************/
int metis_rcode(int sigrval)
{
  switch (sigrval) {
    case 0:
      return METIS_OK;
      break;
    case SIGMEM:
      return METIS_ERROR_MEMORY;
      break;
    default:
      return METIS_ERROR;
      break;
  }
}


