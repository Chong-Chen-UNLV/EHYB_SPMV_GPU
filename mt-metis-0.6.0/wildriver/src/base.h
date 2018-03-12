/**
 * @file base.h
 * @brief Base internal header file.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 *
 */




#ifndef WILDRIVER_BASE_H
#define WILDRIVER_BASE_H




#include "wildriver.h"




/* rename types for internal usage */
#define dim_t wildriver_dim_t
#define ind_t wildriver_ind_t
#define val_t wildriver_val_t




static dim_t const NULL_DIM = (dim_t)-1;
static ind_t const NULL_IND = (ind_t)-1;




#endif
