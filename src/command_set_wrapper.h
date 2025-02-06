#ifndef COMMAND_SET_WRAPPER_H
#define COMMAND_SET_WRAPPER_H

#include "global_basic.h"
#include "command_dist.h"
#include <sys/stat.h>
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <argp.h>
#include <argz.h>
#include <err.h>
#include <errno.h>
#include <math.h>

typedef struct set_opt
{
  int operation;
  bool q2markerdb;
  int p;
  int P;
  int num_remaining_args;
  char ** remaining_args;
  char insketchpath[PATHLEN];
  char pansketchpath[PATHLEN];
  char subsetf[PATHLEN]; //subset infile
  char outdir[PATHLEN];
} set_opt_t ;




int cmd_set(struct argp_state* state);
#endif 
