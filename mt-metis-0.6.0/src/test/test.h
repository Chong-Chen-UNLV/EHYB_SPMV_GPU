/**
 * @file test.h
 * @brief Unit test frame work macros and headers.
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2016
 * @version 1
 * @date 2016-07-18
 */

#ifndef TEST_H
#define TEST_H

#define TESTEQUALS(a,b,fmt) \
  do { \
    if ((a) != (b)) { \
      eprintf("[TEST] " #a " (" fmt ") != " #b " (" fmt ") at %s:%d\n", \
          a,b,__FILE__,__LINE__); \
      return 1; \
    } \
  } while(0)

#define OMPTESTEQUALS(a,b,fmt,r) \
  do { \
    if ((a) != (b)) { \
      eprintf("[TEST]<tid=%d> " #a " (" fmt ") != " #b " (" fmt \
          ") at %s:%d\n",omp_get_thread_num(),a,b,__FILE__,__LINE__); \
      _Pragma("omp atomic") \
      ++(r); \
    } \
  } while(0)

#define TESTLESSTHANOREQUAL(a,b,fmt) \
  do { \
    if ((a) > (b)) { \
      eprintf("[TEST] " #a " (" fmt ") ! <= " #b " (" fmt ") at %s:%d\n", \
          a,b,__FILE__,__LINE__); \
      return 1; \
    } \
  } while(0)

#define TESTGREATERTHANOREQUAL(a,b,fmt) \
  do { \
    if ((a) < (b)) { \
      eprintf("[TEST] " #a " (" fmt ") ! >= " #b " (" fmt ") at %s:%d\n", \
          a,b,__FILE__,__LINE__); \
      return 1; \
    } \
  } while(0)

#define TESTLESSTHAN(a,b,fmt) \
  do { \
    if ((a) >= (b)) { \
      eprintf("[TEST] " #a " (" fmt ") ! < " #b " (" fmt ") at %s:%d\n", \
          a,b,__FILE__,__LINE__); \
      return 1; \
    } \
  } while(0)

#define TESTGREATERTHAN(a,b,fmt) \
  do { \
    if ((a) <= (b)) { \
      eprintf("[TEST] " #a " (" fmt ") ! > " #b " (" fmt ") at %s:%d\n", \
          a,b,__FILE__,__LINE__); \
      return 1; \
    } \
  } while(0)

#define TESTTRUE(a) \
  do { \
    if (!(a)) { \
      eprintf("[TEST] " #a " evaluated to false at %s:%d\n",__FILE__, \
          __LINE__); \
      return 1; \
    } \
  } while (0)

#define OMPTESTTRUE(a,r) \
  do { \
    if (!(a)) { \
      eprintf("[TEST] " #a " evaluated to false at %s:%d\n",__FILE__, \
          __LINE__); \
      _Pragma("omp atomic") \
      ++(r); \
    } \
  } while (0)

int test(void);

#endif
