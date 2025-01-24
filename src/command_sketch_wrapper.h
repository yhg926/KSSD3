#ifndef SKETCH_WRAPPER_H
#define SKETCH_WRAPPER_H

#include "global_basic.h"
#include <stdbool.h>
#include <argp.h>

typedef struct sketch_opt 
{
	// sketch stats
  int hclen; // half context length, 1..16
  int holen; // half outer object length, 0..8
	int iolen; // half outer object length, 0..8
  int drfold; //dimension reduction fold 2^n , 0..32
	int kmerocrs;
	int p ; //threads counts
	bool abundance;
	bool merge_comblco;
	char *fpath;
	char *outdir; // results dir
	char *pipecmd; //pipe command
	int num_remaining_args;
  char **remaining_args;
} sketch_opt_t;

int cmd_sketch(struct argp_state* state);
extern void compute_sketch(sketch_opt_t*, infile_tab_t*);
extern int merge_comblco (sketch_opt_t * sketch_opt_val);

#endif
