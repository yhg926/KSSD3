#ifndef COMMAND_SET
#define COMMAND_SET

#include "command_set_wrapper.h"
//#include "command_operate.h"
typedef struct subset
{
  int taxid;
  char *taxname;
  int *gids;
} subset_t ;

typedef struct compan
{
  int taxn;
  int gn;
  subset_t * tax;
} compan_t;

/*core functions*/
int sketch_union(set_opt_t *set_opt);
int sketch_operate(set_opt_t* set_opt);
int uniq_sketch_union(set_opt_t* set_opt);
int combin_pans(set_opt_t* set_opt);
void print_gnames(set_opt_t* set_opt);
compan_t *organize_taxf(char* taxfile);
int combin_subset_pans(set_opt_t* set_opt);
int grouping_genomes(set_opt_t* set_opt) ; //20230524

#endif 
