#ifndef COMMAND_TRIANGLE
#define COMMAND_TRIANGLE
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

typedef struct triangle_opt
{
	int metric; //0 mashD; 1, aafD; 2. 1-ani
	double c; //minimal distance to enrolled sketches
	int p; //threads
	bool d; //diagnal
  char qrydir[PATHLEN];
  char outf[PATHLEN];
	char gl[PATHLEN]; // genome list with selection code 	
  int num_remaining_args;
  char **remaining_args;

} triangle_opt_t;

int cmd_triangle(struct argp_state* state);
int compute_triangle(triangle_opt_t *);
#endif
