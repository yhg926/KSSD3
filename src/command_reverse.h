#ifndef COMMAND_REVERSE
#define COMMAND_REVERSE
#include "global_basic.h"
#include <stdbool.h>
typedef struct reverse_opt_val
{
	char shufile[PATHLEN];
	char outdir[PATHLEN];
	int p;
	bool byreads;
	int num_remaining_args;
  char **remaining_args;
} reverse_opt_val_t;




int generic_co_reverse2kmer(reverse_opt_val_t *opt_val);
int co_reverse2kmer(reverse_opt_val_t *opt_val);
int co_rvs2kmer_byreads( reverse_opt_val_t *opt_val );
llong core_reverse2unituple(unsigned int kid, int compid, int compbit, int pf_bits, int inner_ctx_bits, int half_outer_ctx_bits, unsigned int *rev_shuf_arr);
#include <argp.h>
int cmd_reverse(struct argp_state* state);



#include <stdint.h> 
static inline void unituple2kstring(uint64_t unituple, int k, char *kstring, char *rc_kstring)
{
	for (int i = k - 1; i >=0; --i)
	{
		kstring[i] = Mapbase[unituple & 0x3u];
		rc_kstring[k - 1 - i] = Mapbase[(~unituple) & 0x3u];
		unituple >>= 2;
	}
	kstring[k] = '\0';
	rc_kstring[k] = '\0';
} // unituple2kstring: convert unituple to kmer string using	

#endif
