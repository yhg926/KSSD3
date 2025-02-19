#ifndef COMMAND_MATRIX
#define COMMAND_MATRIX
#include <argp.h>
#include <math.h>
#include <tgmath.h>
#include "global_basic.h"

// distance macro
#define MIN(X,Y) ( (X) < (Y) ? (X): (Y) )
#define U(X,Y,XnY) ( (X) + (Y) - (XnY))  
#define JCD(X,Y,XnY) ( (double)(XnY)/(U(X,Y,XnY)))
#define CTM(X,Y,XnY) ((double)(XnY)/MIN(X,Y))

#define MASHD (K,X,Y,XnY) (-log( 2*JCD(X,Y,XnY) / (1 + JCD(X,Y,XnY) )) / (K) )
#define AAFD(K,X,Y,XnY) (-log(CTM(X,Y,XnY)) / (K))

typedef struct matrix_opt
{
	int metric; //0 mashD; 1, aafD; 2. 1-ani
	double c; //minimal distance to enrolled sketches
	int p; //threads
	bool d; //diagnal
	bool ani;
	int e;
  char qrydir[PATHLEN];
	char refdir[PATHLEN];
  char outf[PATHLEN];
	char gl[PATHLEN]; // genome list with selection code 	
  int num_remaining_args;
  char **remaining_args;

} matrix_opt_t;

int cmd_matrix(struct argp_state* state);
int compute_triangle(matrix_opt_t *);
int compute_matrix(matrix_opt_t *);
int compute_ani_matrix(matrix_opt_t *matrix_opt);
#endif
