#ifndef COMMAND_OPERATE_H
#define COMMAND_OPERATE_H

#include "command_set.h"

//extern const char sketch_stat[];
//extern const char idx_sketch_suffix[];
void print_lco_gnames(set_opt_t* set_opt);
int lsketch_union(set_opt_t* set_opt);
int lsketch_operate(set_opt_t* set_opt);
int lgrouping_genomes(set_opt_t* set_opt);
#endif
