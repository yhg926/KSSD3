#ifndef SKETCH_WRAPPER_H
#define SKETCH_WRAPPER_H

#include "global_basic.h"
#include <stdbool.h>
#include <argp.h>

typedef struct sketch_opt
{
	// sketch stats
	int hclen;	// half context length, 1..16
	int holen;	// half outer object length, 0..8
	int iolen;	// half outer object length, 0..8
	int drfold; // dimension reduction fold 2^n , 0..32
	int kmerocrs;
	int p; // threads counts
	bool abundance;
	bool conflict; // keep conflict context-object or not 
	bool merge_comblco;
	bool split_mfa;
	bool coden_ctxobj_pattern; 	 
	char index[PATHLEN];
	char *fpath;
	char *outdir;  // results dir
	char *pipecmd; // pipe command
	int num_remaining_args;
	char **remaining_args;
} sketch_opt_t;

int cmd_sketch(struct argp_state *state);
/*
extern void compute_sketch(sketch_opt_t*, infile_tab_t*);
extern void gen_inverted_index4comblco(const char* sketchdir);
extern int merge_comblco (sketch_opt_t * sketch_opt_val);
extern uint32_t get_sketching_id(uint32_t hclen, uint32_t holen,uint32_t iolen,uint32_t drfold,uint32_t FILTER);
*/
#endif
