/**
 * @file test.c
 * @brief Unit test framework.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2016
 * @version 1
 * @date 2016-07-18
 */

#include "test.h"
#include "dlutil.h"

int main(const int argc, const char ** argv) 
{
  dl_timer_t tmr;
  dl_init_timer(&tmr);
  dl_start_timer(&tmr);
  int rv = test();
  double time = dl_poll_timer(&tmr);

  if (rv == 0) {
    printf("PASSED %7.3f s\n",time);
  } else { 
    printf("FAILED %7.3f s\n",time);
  }

  return rv; 
}
