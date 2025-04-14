/*command_operate.c extent the legency kssd set functions to support long sketch (*.comblco) */
#include "global_basic.h"
#include "command_ani.h"
#include "command_operate.h"
#include "../klib/khash.h"
#include <stdint.h>

const char lpan_prefix[]="lpan"; //uint64_t pan
const char luniq_pan_prefix[]="luniq_pan";
// common vars
static size_t file_size; 
static dim_sketch_stat_t lco_stat_readin, lco_stat_pan, lco_stat_origin ;
static struct stat s;
static char outfpath[PATHLEN+20];
static int ret;
extern const char sorted_comb_ctxgid64obj32[];
void print_lco_gnames(set_opt_t* set_opt){
	void *mem_stat = read_from_file( test_get_fullpath(set_opt->insketchpath,sketch_stat) , &file_size);
	memcpy(&lco_stat_readin,mem_stat,sizeof(lco_stat_readin));	char (*tmpname)[PATHLEN] = mem_stat+sizeof(dim_sketch_stat_t);
	uint64_t *mem_index = (uint64_t *)read_from_file( test_get_fullpath(set_opt->insketchpath,idx_sketch_suffix) , &file_size);
	for(int i = 0 ; i<lco_stat_readin.infile_num ; i++) printf("%lu\t%s\n",mem_index[i+1]-mem_index[i],tmpname[i]);
	free_all(mem_stat,mem_index,NULL)	;
}

void show_content(set_opt_t* set_opt){
	FILE *fp; uint64_t kmer; ctxgidobj_t ctxgidobj;
	if(set_opt->show == 1) {
		uint64_t *sketch_index = read_from_file( test_get_fullpath(set_opt->insketchpath,idx_sketch_suffix) , &file_size);
		int infile_num = file_size/sizeof(uint64_t) - 1;	
		if((fp=fopen(test_get_fullpath(set_opt->insketchpath,combined_sketch_suffix),"rb")) == NULL) 
			err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'", __func__, set_opt->insketchpath,combined_sketch_suffix);
		for(int i = 0 ; i< infile_num;i++){
			for(uint64_t j = sketch_index[i]; j < sketch_index[i+1]; j++ ){
			fread(&kmer,sizeof(kmer),1,fp); printf("%d\t%lx\n",i,kmer);
			}
		}
		free(sketch_index);	
	}
	else if (set_opt->show == 2){
		if((fp=fopen(test_get_fullpath(set_opt->insketchpath,sorted_comb_ctxgid64obj32),"rb")) == NULL)
            err(EXIT_FAILURE, "%s(): Failed to open file '%s/%s'", __func__, set_opt->insketchpath,sorted_comb_ctxgid64obj32);
		while(!feof(fp)){
			fread(&ctxgidobj,sizeof(ctxgidobj),1,fp);  printf("%lx\t%x\n",ctxgidobj.ctxgid,ctxgidobj.obj);
		}
	}
	else err(EXIT_FAILURE, "%s(): only set_opt.show values 1 and 2 are supported, set_opt.show =%d", __func__, set_opt->show);
	fclose(fp);
}

KHASH_MAP_INIT_INT64(kmer_hash, int)

int lsketch_union(set_opt_t* set_opt){ // for both union and uniq union

	void *mem_stat = read_from_file( test_get_fullpath(set_opt->insketchpath,sketch_stat) , &file_size);
	memcpy(&lco_stat_readin,mem_stat,sizeof(lco_stat_readin));
  	if(lco_stat_readin.infile_num == 1){ // no need create 
    	printf("only 1 sketch, use %s as pan-sketch?(Y/N)\n",set_opt->insketchpath);
    	char inpbuff; scanf(" %c", &inpbuff);
    	if ( (inpbuff == 'Y') || (inpbuff == 'y') ) {
      		chdir(set_opt->insketchpath);
      		if(rename(combined_sketch_suffix,lpan_prefix) !=0)  err(errno,"lsketch_union()");
      		printf("the union directory: %s created successfully\n", set_opt->insketchpath) ;
      	return 1;
    	}
  	}
 	write_to_file(test_create_fullpath(set_opt->outdir,sketch_stat), mem_stat,file_size);
//union operation
	uint64_t *mem_comblco = (uint64_t *)read_from_file( test_get_fullpath(set_opt->insketchpath,combined_sketch_suffix) , &file_size);
	uint32_t in_kmer_ct = file_size/sizeof(mem_comblco[0]);
	khash_t(kmer_hash) *h = kh_init(kmer_hash); //kh_init_size(kmer_hash, (float)in_kmer_ct * 1.66 );
	for(uint32_t i = 0; i< in_kmer_ct; i++){
		khint_t key = kh_put(kmer_hash, h, mem_comblco[i], &ret);
		if(ret) kh_value(h, key) = 1;
		else kh_value(h, key)++;
	}
	//uint32_t pan_kmer_ct = 0;
	Vector mem_pan; vector_init(&mem_pan, sizeof(uint64_t));
	if(set_opt->operation == 2){ // -u: normal union mode
		for (khint_t k = kh_begin(h); k != kh_end(h); ++k) {
  			if (kh_exist(h, k))  vector_push(&mem_pan, &kh_key(h, k)); //mem_pan[pan_kmer_ct++] = kh_key(h, k);   
  		}
		write_to_file(test_create_fullpath(set_opt->outdir,lpan_prefix),mem_pan.data, mem_pan.element_size * mem_pan.size);
	}
	else if(set_opt->operation == 3){ // -q: uniq union mode	
		if(!set_opt->q2markerdb){
			for (khint_t k = kh_begin(h); k != kh_end(h); ++k) {
                if (kh_exist(h, k) && kh_value(h, k) == 1)  vector_push(&mem_pan, &kh_key(h, k));
            }
            write_to_file(test_create_fullpath(set_opt->outdir,luniq_pan_prefix),mem_pan.data, mem_pan.element_size * mem_pan.size);
		}
		else{ // genereate markerdb directly instead of uniq union 
			uint64_t *fco_pos = (uint64_t *)read_from_file( test_get_fullpath(set_opt->insketchpath,idx_sketch_suffix),&file_size);
    		uint64_t *post_fco_pos = calloc((lco_stat_readin.infile_num + 1),sizeof(uint64_t));
		    for (uint32_t i = 0; i <lco_stat_readin.infile_num;i++){
        		for(uint64_t n = fco_pos[i]; n < fco_pos[i+1] ; n++){
					khint_t k = kh_get(kmer_hash, h, mem_comblco[n]) ;
            		if(kh_value(h, k) == 1) vector_push(&mem_pan,&mem_comblco[n]);           
        		}
        		post_fco_pos[i+1] = mem_pan.size;
    		}
    		write_to_file(format_string("%s/%s",set_opt->outdir,combined_sketch_suffix), mem_pan.data, mem_pan.element_size * mem_pan.size);
		    write_to_file(format_string("%s/%s",set_opt->outdir,idx_sketch_suffix),post_fco_pos,(lco_stat_readin.infile_num+1) * sizeof(post_fco_pos[0]));
			free_all(fco_pos,post_fco_pos,NULL);
		}
	}
	else err(EINVAL,"operation value %d neither 2 (-u: union) nor 3 (-q :uniq uion )",set_opt->operation);
	kh_destroy(kmer_hash, h);  vector_free(&mem_pan);
	free_all(mem_stat,mem_comblco,NULL);	
  return 1;
}

KHASH_SET_INIT_INT64(kmer_set)
int lsketch_operate(set_opt_t* set_opt){
clock_t start_time = clock();
	void *mem_stat_pan = read_from_file( test_get_fullpath(set_opt->pansketchpath,sketch_stat) , &file_size);
	memcpy(&lco_stat_pan,mem_stat_pan,sizeof(lco_stat_pan));
	void *mem_stat_lco = read_from_file( test_get_fullpath(set_opt->insketchpath,sketch_stat), &file_size);
	memcpy(&lco_stat_origin, mem_stat_lco,sizeof(lco_stat_origin));
	if (lco_stat_pan.hash_id != lco_stat_origin.hash_id) err(EXIT_FAILURE,"%s(): %s sketcing id %u != %s id %u",__func__,set_opt->pansketchpath, lco_stat_origin.hash_id, set_opt->insketchpath, lco_stat_pan.hash_id);	
	if(lco_stat_origin.koc) printf("%s() Warning: k-mer abundances are dropped in this sketch operation\n ",__func__);
//copy sketch stat file to result sketch
	write_to_file(test_create_fullpath(set_opt->outdir,sketch_stat),mem_stat_lco,file_size);

	char *lco_fpath; 
	if( (lco_fpath = test_get_fullpath(set_opt->pansketchpath,lpan_prefix)) == NULL && (lco_fpath = test_get_fullpath(set_opt->pansketchpath,luniq_pan_prefix)) == NULL  )
		err(EXIT_FAILURE,"%s():cannot find %s or %s under %s ",__func__,lpan_prefix,luniq_pan_prefix, set_opt->pansketchpath);
	uint64_t *mem_pan = (uint64_t *)read_from_file(lco_fpath,&file_size);			
  khash_t(kmer_set) *h = kh_init(kmer_set); 
  uint32_t kmer_ct = file_size/sizeof(uint64_t);
	for(int i = 0 ; i< kmer_ct; i++) kh_put(kmer_set,h,mem_pan[i],&ret);
  //read comblco.index to mem	
	uint64_t *fco_pos = (uint64_t *)read_from_file( test_get_fullpath(set_opt->insketchpath,idx_sketch_suffix),&file_size);
	// post operation index	calloc	 
 	uint64_t *post_fco_pos = calloc((lco_stat_origin.infile_num + 1),sizeof(uint64_t));	
	//read comblco to mem
	uint64_t *tmp_comblco_mem = (uint64_t *)read_from_file( test_get_fullpath(set_opt->insketchpath,combined_sketch_suffix),&file_size);
	uint64_t *post_comblco_mem = (uint64_t *) malloc(file_size); uint32_t post_kmer_ct = 0;
   
	// sketch operation
	for (uint32_t i = 0; i <lco_stat_origin.infile_num;i++){
		for(uint64_t n = fco_pos[i]; n < fco_pos[i+1] ; n++){
			 //make sure set_opt->operation == 0 if subtract, == 1 if intersect
			if(set_opt->operation == (kh_get(kmer_set, h, tmp_comblco_mem[n]) != kh_end(h)) )
				post_comblco_mem[post_kmer_ct++] = tmp_comblco_mem[n];				
		}
		post_fco_pos[i+1] = post_kmer_ct;
	}
	kh_destroy(kmer_set,h);
 //write to result comblco
  sprintf(outfpath,"%s/%s",set_opt->outdir,combined_sketch_suffix);
	write_to_file(outfpath,post_comblco_mem,post_kmer_ct * sizeof(post_comblco_mem[0]));
//write index
	sprintf(outfpath,"%s/%s",set_opt->outdir,idx_sketch_suffix);
	write_to_file(outfpath,post_fco_pos,(lco_stat_origin.infile_num+1) * sizeof(post_fco_pos[0]));

	free_all(mem_stat_pan, mem_stat_lco,lco_fpath,mem_pan,fco_pos, post_fco_pos, tmp_comblco_mem,post_comblco_mem,NULL);
	return 1;
}



int lgrouping_genomes(set_opt_t* set_opt){

	void *mem_stat = read_from_file( test_get_fullpath(set_opt->insketchpath,sketch_stat) , &file_size);
	memcpy(&lco_stat_readin,mem_stat,sizeof(lco_stat_readin));
		
	compan_t *subset = organize_taxf(set_opt->subsetf);
	if(lco_stat_readin.infile_num != subset->gn)  err(EXIT_FAILURE,"%s(): %s's genome number %d != %s's line number %d",__func__,set_opt->insketchpath,lco_stat_readin.infile_num,set_opt->subsetf,subset->gn);
 
//read index and comblco
	uint64_t *tmp_idx = (uint64_t *)read_from_file( test_get_fullpath(set_opt->insketchpath,idx_sketch_suffix), &file_size);
	uint64_t *mem_comblco = (uint64_t *)read_from_file( test_get_fullpath(set_opt->insketchpath,combined_sketch_suffix), &file_size);

	if(tmp_idx[subset->gn]*sizeof(uint64_t) != file_size) err(EXIT_FAILURE,"%s(): %s last(%u) index(%lu) * sizeof(uint64_t) != %s file size (%lu) ", \
			 __func__,idx_sketch_suffix,subset->gn,tmp_idx[subset->gn], combined_sketch_suffix ,file_size);
// out index and comblco
	uint64_t *grouped_comblco = (uint64_t *)malloc(file_size);  	
	uint64_t *out_idx  = calloc( (subset->taxn+1), sizeof(uint64_t));
  int outfn = 0; uint64_t grouped_kmer_ct = 0;

	for(int t = 0; t < subset->taxn; t++){
		if(subset->tax[t].taxid == 0) continue;// ignore taxid 0		
		outfn++; out_idx[outfn] = out_idx[outfn-1];
		khash_t(kmer_set) *h = kh_init(kmer_set); 
		for(int n = 1; n <= subset->tax[t].gids[0];n++){
			int gid = subset->tax[t].gids[n] ;		
			for(uint64_t i = tmp_idx[gid]; i < tmp_idx[gid+1]; i++)
      	kh_put(kmer_set, h, mem_comblco[i],&ret);		
		}
		//write into comblco for tax t
		for (khint_t k = kh_begin(h); k != kh_end(h); ++k) {
      if (kh_exist(h, k)) {
				grouped_comblco[grouped_kmer_ct++] = kh_key(h, k);
				out_idx[outfn]++;	
			}
    }	
		kh_destroy(kmer_set,h);		
	}
	//write grouped kmer and index to result
	write_to_file(test_create_fullpath(set_opt->outdir,combined_sketch_suffix), grouped_comblco, grouped_kmer_ct*sizeof(grouped_comblco[0]));
  write_to_file(test_create_fullpath(set_opt->outdir,idx_sketch_suffix), out_idx, (outfn+1)*sizeof(out_idx[0]));
	//write stat file 
	lco_stat_readin.infile_num = outfn; lco_stat_readin.koc = 0;
	char (*tmpfname)[PATHLEN] = malloc(outfn * PATHLEN); int idx = 0;
  for(int t=0;t<subset->taxn;t++){
    if(subset->tax[t].taxid != 0) {
      if(subset->tax[t].taxname != NULL)
        sprintf(tmpfname[idx],"%d_%s",subset->tax[t].taxid,subset->tax[t].taxname);
      else sprintf(tmpfname[idx],"%d",subset->tax[t].taxid);
      idx++;
    }
    free_all(subset->tax[t].gids, subset->tax[t].taxname,NULL);
  }
	concat_and_write_to_file(test_create_fullpath(set_opt->outdir,sketch_stat),&lco_stat_readin,sizeof(lco_stat_readin),tmpfname,PATHLEN*outfn);
	
	free_all( mem_stat, mem_comblco,grouped_comblco,tmp_idx,out_idx,subset->tax,subset,tmpfname,NULL );
 	return outfn;
}

































