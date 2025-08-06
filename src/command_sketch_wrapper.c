#include "command_sketch_wrapper.h"
#include "sketch_rearrange.h"
#include "global_basic.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <argp.h>
#include <argz.h>
#include <sys/stat.h>
#include <dirent.h>
#include <err.h>
#include <errno.h>
#include <unistd.h>
#include <stdbool.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

/*** argp wrapper ***/
struct arg_sketch
{
  struct arg_global *global;

  char *name;
};

static struct argp_option opt_sketch[] =
    {
        {"ctxlen", 'C', "<INT>", 0, "Half context length. [11]\v", 1},
        {"outerobjlen", 'O', "<INT>", 0, "Half outer object length. [0]\v", 2},
        {"innerobjlen", 'I', "<INT>", 0, "Inner object length. [0]\v", 3},
        {"DimRdcFold", 'f', "<INT>", 0, "K-mer space downsampling rate 1/2^f. [12]\v", 4},
        {"LstKmerOcrs", 'n', "INT", 0, "Specify the Least Kmer occurence in fastq file. [1]\v", 5},
        {"list", 'l', "file", 0, "a file contain paths for all query sequences\v", 6},
        {"outdir", 'o', "<path>", 0, "folder path for results files.\v", 7},
        {"abundance", 'A', 0, 0, "abundance estimate mode.\v", 8},
        {"use_coden_ctxobj", 'T', 0, 0, "use coden context object structure pattern [0].\v", 8},
        {"threads", 'p', "<INT>", 0, "Threads number to use. [1]\v", 9},
        {"index", 'i', "<FILE>", 0, "Inverted indexing of combined sketch.\v", 10},
        {"conflict", 666, 0, 0, "keep conflict context-objet(for raw reads sketching).\v", 10},
        {"merge", 777, 0, 0, "merge sketches, not for genome sketching.\v", 10},
        {"splitmfa", 888, 0, 0, "treat mfa file as many genomes.\v", 11},
        {0}  
    };
static char doc_sketch[] =
    ""
    ""
    "\v"
    "";

sketch_opt_t sketch_opt = {
    .hclen = 11, //
    .holen = 0,
    .iolen = 0,
    .drfold = 12,
    .kmerocrs = 1,
    .p = 1,         // threads num: p
    .abundance = 0, // no abundance
    .conflict = 0, // no conflict context-objet.
    .merge_comblco = 0,
    .split_mfa = 0,
    .coden_ctxobj_pattern = false, // no coden ctxobj pattern
                                   //  .fpath[0] ='\0',
    .outdir = "./",
    .index[0] = '\0',
    //	.pipecmd[0] = '\0', // no pipe command
    .num_remaining_args = 0, // int num_remaining_args; no option arguments num.
    .remaining_args = NULL   // char **remaining_args; no option arguments array.
};

static error_t parse_sketch(int key, char *arg, struct argp_state *state)
{
  struct arg_sketch *sketch = state->input;
  switch (key)
  {
  case 'C':
  {
    int val = atoi(arg);
    if (val < 4 || val > 16)
    {
      printf("-C requires a integer within ranges 4..16 \n");
      exit(1);
    }
    sketch_opt.hclen = val;
    break;
  }
  case 'O':
  {
    int val = atoi(arg);
    if (val < 0 || val > 8)
    {
      printf("-O requires a integer within ranges 0..8 \n");
      exit(1);
    }
    sketch_opt.holen = val;
    break;
  }
  case 'I':
  {
    int val = atoi(arg);
    if (val < 0 || val > 8)
    {
      printf("-O requires a integer within ranges 0..8 \n");
      exit(1);
    }
    sketch_opt.iolen = val;
    break;
  }
  case 'f':
  {
    int val = atoi(arg);
    if (val < 0 || val > 24)
    {
      printf("-f requires a integer within ranges 0..24 \n");
      exit(1);
    }
    sketch_opt.drfold = val;
    break;
  }
  case 'n':
  {
    int val = atoi(arg);
    if (val < 1 || val > 65536)
      err(errno, "%s():-n requires a integer within ranges 1..65536", __func__);

    sketch_opt.kmerocrs = val;
    break;
  }
  case 'i':
  {
    strcpy(sketch_opt.index, arg);
    break;
  }
  case 'p':
  {
#ifdef _OPENMP
    sketch_opt.p = atoi(arg);
#else
    warnx("This version of kssd was built without OpenMP and "
          "thus does not support multi threading. Ignoring -p %d",
          atoi(arg));
#endif
    break;
  }
  case 'A':
  {
    sketch_opt.abundance = 1;
    break;
  }
  case 'T':
  {
    sketch_opt.coden_ctxobj_pattern = 1;
    break;
  }
  case 'P':
  {
    sketch_opt.pipecmd = malloc(strlen(arg));
    strcpy(sketch_opt.pipecmd, arg);
    break;
  }
  case 'l':
  {
    sketch_opt.fpath = malloc(strlen(arg));
    strcpy(sketch_opt.fpath, arg);
    break;
  }
  case 'o':
  {
    sketch_opt.outdir = malloc(strlen(arg) + 10);
    strcpy(sketch_opt.outdir, arg);
    break;
  }
  case 666:
  {
    sketch_opt.conflict = 1; // keep conflict context-objet(for raw reads sketching).
    break;
  }
  case 777:
  {
    sketch_opt.merge_comblco = 1;
    break;
  }
  case 888:
  {
    sketch_opt.split_mfa = 1;
    break;
  }
  case ARGP_KEY_ARGS:
  {
    sketch_opt.num_remaining_args = state->argc - state->next;
    sketch_opt.remaining_args = state->argv + state->next;
    break;
  }
  case ARGP_KEY_END:
  {
    int klen = sketch_opt.iolen + 2 * (sketch_opt.holen + sketch_opt.hclen);
    if (klen > 32)
    {
      printf("\nError: k-mer length %d should smaller than 32 \n\n", klen);
      exit(1);
    }
    if (sketch_opt.index[0] == '\0' && sketch_opt.fpath == NULL && sketch_opt.num_remaining_args == 0)
    {
      printf("\nError: missing input sequences file \n\n");
      argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
      argp_usage(state);
    }
    break;
  }
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp_sketch =
    {
        opt_sketch,
        parse_sketch,
        0, //  "[arguments ...]",
        doc_sketch};

infile_tab_t *sketch_organize_infiles(sketch_opt_t *sketch_opt_val)
{
  int fmt_ck;
  if (sketch_opt_val->pipecmd == NULL)
    fmt_ck = 1; // need check format-- normal mode
  else
    fmt_ck = 0;
  // do it if fpath is not ""
  if (sketch_opt_val->fpath != NULL)
  {
    return organize_infile_list(sketch_opt_val->fpath, fmt_ck);
  }
  else if (sketch_opt_val->num_remaining_args > 0)
  {
    return organize_infile_frm_arg(sketch_opt_val->num_remaining_args, sketch_opt_val->remaining_args, fmt_ck);
  }
  else
  {
    perror("please specify the (meta)genome files");
    return NULL;
  }
};

extern uint32_t FILTER; // control dimensionality reducation level
extern uint32_t hash_id;
extern dim_sketch_stat_t comblco_stat_one;
extern void compute_sketch(sketch_opt_t *, infile_tab_t *);
extern void gen_inverted_index4comblco(const char *sketchdir);
extern int merge_comblco(sketch_opt_t *sketch_opt_val);
extern uint32_t get_sketching_id(uint32_t hclen, uint32_t holen, uint32_t iolen, uint32_t drfold, uint32_t FILTER);

int cmd_sketch(struct argp_state *state)
{
  struct arg_sketch sketch = {
      0,
  };
  int argc = state->argc - state->next + 1;
  char **argv = &state->argv[state->next - 1];
  sketch.global = state->input;
  argp_parse(&argp_sketch, argc, argv, ARGP_IN_ORDER, &argc, &sketch);
  state->next += argc - 1;

  if (sketch_opt.merge_comblco)
  {
    int merge_count = merge_comblco(&sketch_opt);
  }
  else if (sketch_opt.index[0] != '\0')
  {
    gen_inverted_index4comblco(sketch_opt.index);
  }
  else
  {
    infile_tab_t *infile_stat = sketch_organize_infiles(&sketch_opt);
    if (infile_stat->infile_num)
    {
      FILTER = UINT32_MAX >> sketch_opt.drfold;
      /* conditionally initilize some comblco_stat_one members*/
      if(sketch_opt.coden_ctxobj_pattern){
        hash_id = get_sketching_id(NUM_CODENS,0,0,sketch_opt.drfold, FILTER);
        klen = 3 * NUM_CODENS + 1; // klen is 3*NUM_CODENS+1
        comblco_stat_one.coden_len = NUM_CODENS; // set coden_len
        comblco_stat_one.hclen = 0;
        comblco_stat_one.holen = 0;
      }else{
        hash_id =  get_sketching_id(sketch_opt.hclen, sketch_opt.holen, sketch_opt.iolen, sketch_opt.drfold, FILTER);
        klen = 2 * (sketch_opt.hclen + sketch_opt.holen) + sketch_opt.iolen;
        comblco_stat_one.coden_len = 0; // no coden ctxobj pattern
        comblco_stat_one.hclen = sketch_opt.hclen;
        comblco_stat_one.holen = sketch_opt.holen;
      }
      if (klen > 32 || FILTER < 256)
        err(EINVAL, "%s(): klen (%d) or FILTER (%u) is out of range (klen <=32 and FILTER: 256..0xffffffff)", __func__, klen, FILTER);

      printf("Sketching method hashid = %u\tctxobj_coden_len=%u\tklen=%u\tFILTER=%u\thclen=%d\n", hash_id, comblco_stat_one.coden_len, klen, FILTER, comblco_stat_one.hclen);
      { /*initilize the rest comblco_stat_one member*/
        comblco_stat_one.hash_id = hash_id;
        comblco_stat_one.koc = sketch_opt.abundance;
        comblco_stat_one.klen = klen;     
        comblco_stat_one.drfold = sketch_opt.drfold;
        comblco_stat_one.infile_num = infile_stat->infile_num;
      }
      // ensure initalize global masks
      const_comask_init(&comblco_stat_one);
      mkdir_p(sketch_opt.outdir);
      // initialize k-mer rearrange method: reorder_unituple_by_coden_pattern64() or uint64_kmer2ctxobj()
      set_uint64kmer2generic_ctxobj(sketch_opt.coden_ctxobj_pattern);
      compute_sketch(&sketch_opt, infile_stat);
    }
    else
    {
      printf("not valid fas/fastq files!\n");
    }
    free(infile_stat->organized_infile_tab);
  }
  return 1;
};
