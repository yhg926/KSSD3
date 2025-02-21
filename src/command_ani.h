#ifndef COMMAND_ANI
#define COMMAND_ANI
#include <argp.h>
#include <math.h>
#include <tgmath.h>
#include "global_basic.h"

typedef struct ani_opt
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

} ani_opt_t;


typedef uint32_t obj_t ; // optional uint64_t
typedef struct id_obj
{
	uint32_t gid;
	obj_t obj;
} id_obj_t;

int cmd_ani(struct argp_state* state);
int compute_ani(ani_opt_t *ani_opt);
void gen_inverted_index4comblco(const char *refdir);
//int ani_triangle(ani_opt_t *);
//int ani_matrix(ani_opt_t *);
#endif
