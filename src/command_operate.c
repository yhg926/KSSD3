/*command_operate.c extent the legency kssd set functions to support long sketch (*.comblco) */
#include "global_basic.h"
#include "command_operate.h"
#include "klib/khash.h"
#include <stdint.h>



void print_lco_gnames(set_opt_t* set_opt){

  const char* lco_stat_fpath  = test_get_fullpath(set_opt->insketchpath,sketch_stat);
  if(lco_stat_fpath == NULL ) err(errno,"print_lco_gnames(): cannot find %s under %s ",sketch_stat,set_opt->insketchpath);
	const char* lco_index_fpath = test_get_fullpath(set_opt->insketchpath,idx_sketch_suffix);
	if(lco_index_fpath == NULL ) err(errno,"print_lco_gnames(): cannot find %s under %s ",idx_sketch_suffix,set_opt->insketchpath);

  FILE *lco_stat_fp, *lco_index_fp;
  if( ( lco_stat_fp = fopen(lco_stat_fpath,"rb")) == NULL ) err(errno,"print_lco_gnames():%s",lco_stat_fpath);
	if( ( lco_index_fp = fopen(lco_index_fpath,"rb")) == NULL ) err(errno,"print_lco_gnames():%s",lco_index_fpath);

  dim_sketch_stat_t lco_stat_readin;
  fread( &lco_stat_readin, sizeof(lco_stat_readin),1,lco_stat_fp);
  char (*tmpname)[PATHLEN] = malloc( PATHLEN * lco_stat_readin.infile_num );
  fread(tmpname,PATHLEN,lco_stat_readin.infile_num,lco_stat_fp);
	
	uint64_t *tmp_ctx_ct = malloc(sizeof(uint64_t) * (lco_stat_readin.infile_num + 1));
	fread(tmp_ctx_ct,sizeof(uint64_t),(lco_stat_readin.infile_num+1),lco_index_fp) ;

  for(int i=0; i<lco_stat_readin.infile_num; i++ ){
    printf("%lu\t%s\n",tmp_ctx_ct[i+1]-tmp_ctx_ct[i],tmpname[i]);
  }
  free(tmp_ctx_ct);
  free(tmpname);
	fclose(lco_stat_fp);
	fclose(lco_index_fp);
}

const char lpan_prefix[]="lpan"; //uint64_t pan
const char luniq_pan_prefix[]="luniq_pan";

KHASH_MAP_INIT_INT64(kmer_hash, int)
int lsketch_union(set_opt_t* set_opt){ // for both union and uniq union
  const char* lco_stat_fpath = test_get_fullpath(set_opt->insketchpath,sketch_stat);
  if(lco_stat_fpath == NULL ) err(errno,"cannot find %s under %s ",sketch_stat,set_opt->insketchpath);

  FILE *lco_stat_fp;
  if( ( lco_stat_fp = fopen(lco_stat_fpath,"rb")) == NULL ) err(errno,"lsketch_union():%s",lco_stat_fpath);
  dim_sketch_stat_t lco_stat_readin;
  fread( &lco_stat_readin, sizeof(lco_stat_readin),1,lco_stat_fp );

  if(lco_stat_readin.infile_num == 1){ // no need create
    char inpbuff;
    printf("only 1 sketch, use %s as pan-sketch?(Y/N)\n",set_opt->insketchpath);
    scanf(" %c", &inpbuff);

    if ( (inpbuff == 'Y') || (inpbuff == 'y') ) {
      chdir(set_opt->insketchpath);
      if(rename(combined_sketch_suffix,lpan_prefix) !=0)  err(errno,"lsketch_union()");
      printf("the union directory: %s created successfully\n", set_opt->insketchpath) ;
      return 1;
    }
  }
  mkdir_p(set_opt->outdir);
  char outpath[PATHLEN+20];
  sprintf(outpath,"%s/%s",set_opt->outdir,sketch_stat);
  FILE *lco_stat_fp2;
  if( ( lco_stat_fp2 = fopen(outpath,"wb")) == NULL ) err(errno,"lsketch_union():%s",outpath);
  fwrite( &lco_stat_readin,sizeof(lco_stat_readin),1, lco_stat_fp2 );
  fclose(lco_stat_fp);
  fclose(lco_stat_fp2);

	sprintf(outpath,"%s/%s",set_opt->insketchpath,combined_sketch_suffix);		
	FILE *lsketch_fp;
	if( ( lsketch_fp = fopen(outpath,"rb")) == NULL ) err(errno,"lsketch_union():%s",outpath);

	int ret;uint64_t tmpcbdlco;
  khash_t(kmer_hash) *h = kh_init(kmer_hash);
	while (fread(&tmpcbdlco, sizeof(tmpcbdlco),1,lsketch_fp  ) > 0 ){
		khint_t key = kh_put(kmer_hash, h,tmpcbdlco,&ret); 	
		if(ret) kh_value(h, key) = 1;
		else kh_value(h, key)++;
	}
	fclose(lsketch_fp);

	if(set_opt->operation == 2){ // -u: normal union mode
	  sprintf(outpath,"%s/%s",set_opt->outdir,lpan_prefix);
  	if( ( lsketch_fp = fopen(outpath,"wb")) == NULL ) err(errno,"lsketch_union():%s",outpath);
		for (khint_t k = kh_begin(h); k != kh_end(h); ++k) {
  		if (kh_exist(h, k))  fwrite(&kh_key(h, k),sizeof(uint64_t),1,lsketch_fp );    
  	}
		fclose(lsketch_fp);
	}
	else if(set_opt->operation == 3){ // -q: uniq union mode
		sprintf(outpath,"%s/%s",set_opt->outdir,luniq_pan_prefix);
		if( ( lsketch_fp = fopen(outpath,"wb")) == NULL ) err(errno,"lsketch_union():%s",outpath);
	  for (khint_t k = kh_begin(h); k != kh_end(h); ++k) {
      if (kh_exist(h, k) && kh_value(h, k) == 1)  fwrite(&kh_key(h, k),sizeof(uint64_t),1,lsketch_fp );
    }
	  fclose(lsketch_fp);
	}
	else err(EINVAL,"operation value %d neither 2 (-u: union) nor 3 (-q :uniq uion )",set_opt->operation);

	kh_destroy(kmer_hash, h);
		
  return 1;
}

KHASH_SET_INIT_INT64(kmer_set)
int lsketch_operate(set_opt_t* set_opt){
	dim_sketch_stat_t lco_stat_pan, lco_stat_origin ; 
	struct stat s;
	char outfpath[PATHLEN+20];
	char* lco_stat_fpath = test_get_fullpath(set_opt->pansketchpath,sketch_stat);
	if(lco_stat_fpath == NULL ) err(errno,"lsketch_operate():cannot find %s under %s ",sketch_stat,set_opt->pansketchpath);
	FILE *lco_stat_fp;
	if( ( lco_stat_fp = fopen(lco_stat_fpath,"rb")) == NULL ) err(errno,"lsketch_operate():%s",lco_stat_fpath);
	fread( &lco_stat_pan, sizeof(lco_stat_pan),1,lco_stat_fp );
	fclose(lco_stat_fp); free(lco_stat_fpath);

	lco_stat_fpath = test_get_fullpath(set_opt->insketchpath,sketch_stat);
	if(lco_stat_fpath == NULL ) err(errno,"lsketch_operate(): cannot find %s under %s ",sketch_stat,set_opt->insketchpath);
  if( ( lco_stat_fp = fopen(lco_stat_fpath,"rb")) == NULL ) err(errno,"lsketch_operate():%s",lco_stat_fpath);
	stat(lco_stat_fpath, &s);
	void *tmp_stat_mem = malloc(s.st_size);
	fread(tmp_stat_mem,s.st_size,1,lco_stat_fp);
	fclose(lco_stat_fp); free(lco_stat_fpath);

	lco_stat_origin = *(dim_sketch_stat_t*)tmp_stat_mem;
	if (lco_stat_pan.hash_id != lco_stat_origin.hash_id) err(errno,"sketcing id not match(%d Vs. %d)",lco_stat_origin.hash_id,lco_stat_pan.hash_id);	
//copy sketch stat file to result sketch
	mkdir_p(set_opt->outdir) ;
	sprintf(outfpath,"%s/%s",set_opt->outdir,sketch_stat);
	if( ( lco_stat_fp = fopen(outfpath,"wb")) == NULL ) err(errno,"lsketch_operate():%s",outfpath);
	fwrite(tmp_stat_mem,s.st_size,1,lco_stat_fp);
	fclose(lco_stat_fp); 
	free(tmp_stat_mem);
	
	if( (lco_stat_fpath = test_get_fullpath(set_opt->pansketchpath,lpan_prefix)) == NULL
		&& (lco_stat_fpath = test_get_fullpath(set_opt->pansketchpath,luniq_pan_prefix)) == NULL  )
		err(errno,"lsketch_operate():cannot find %s or %s under %s ",lpan_prefix,luniq_pan_prefix, set_opt->pansketchpath);
	
	if( ( lco_stat_fp = fopen(lco_stat_fpath,"rb")) == NULL ) err(errno,"lsketch_operate():%s",lco_stat_fpath);
	
  khash_t(kmer_set) *h = kh_init(kmer_set);
	int ret; uint64_t tmpcbdlco;
  while (fread(&tmpcbdlco, sizeof(tmpcbdlco),1,lco_stat_fp ) > 0 )
    kh_put(kmer_set, h,tmpcbdlco,&ret);
  fclose(lco_stat_fp);free(lco_stat_fpath);

  uint64_t *fco_pos = malloc(sizeof(size_t) * (lco_stat_origin.infile_num + 1) );
 	uint64_t *post_fco_pos = calloc((lco_stat_origin.infile_num + 1),sizeof(size_t));
	//read comblco.index to mem
	lco_stat_fpath = test_get_fullpath(set_opt->insketchpath,idx_sketch_suffix);
	if( ( lco_stat_fp = fopen(lco_stat_fpath,"rb")) == NULL ) err(errno,"lsketch_operate():%s",lco_stat_fpath);
	fread(fco_pos,sizeof(uint64_t),lco_stat_origin.infile_num + 1 ,lco_stat_fp);
	fclose(lco_stat_fp);free(lco_stat_fpath);
	//read comblco to mem
	lco_stat_fpath = test_get_fullpath(set_opt->insketchpath,combined_sketch_suffix);
	if( ( lco_stat_fp = fopen(lco_stat_fpath,"rb")) == NULL ) err(errno,"lsketch_operate():%s",lco_stat_fpath);
  uint64_t *tmp_comblco_mem = malloc(fco_pos[lco_stat_origin.infile_num]*sizeof(uint64_t));
	fread(tmp_comblco_mem,sizeof(uint64_t), fco_pos[lco_stat_origin.infile_num],lco_stat_fp);	
	fclose(lco_stat_fp);free(lco_stat_fpath);

	//write to result comblco
	sprintf(outfpath,"%s/%s",set_opt->outdir,combined_sketch_suffix);
	if( ( lco_stat_fp = fopen(outfpath,"wb")) == NULL ) err(errno,"lsketch_operate():%s",outfpath);
	for (uint32_t i = 0; i <lco_stat_origin.infile_num;i++){
		post_fco_pos[i+1] = post_fco_pos[i];
		for(uint32_t n = 0; n < fco_pos[i+1] - fco_pos[i]; n++){
			 //make sure set_opt->operation == 0 if subtract, == 1 if intersect
			if(set_opt->operation == (kh_get(kmer_set, h, tmp_comblco_mem[fco_pos[i]+n]) != kh_end(h)) ){
				fwrite(tmp_comblco_mem+fco_pos[i]+n,sizeof(uint32_t),1,lco_stat_fp);
				post_fco_pos[i+1]++;
			}
		}
	}
	kh_destroy(kmer_set,h);
	fclose(lco_stat_fp);
	free(fco_pos);free(tmp_comblco_mem);

//write index
	sprintf(outfpath,"%s/%s",set_opt->outdir,idx_sketch_suffix);
  if( ( lco_stat_fp = fopen(outfpath,"wb")) == NULL ) err(errno,"lsketch_operate():%s",outfpath);
	fwrite(post_fco_pos,sizeof(post_fco_pos[0]),lco_stat_origin.infile_num + 1,lco_stat_fp);
	free(post_fco_pos);fclose(lco_stat_fp);
	return 1;
}



int lgrouping_genomes(set_opt_t* set_opt){

	compan_t *subset = organize_taxf(set_opt->subsetf);
	char* lco_stat_fpath = test_get_fullpath(set_opt->insketchpath,sketch_stat);
	FILE *tmpfh;
  if( ( tmpfh = fopen(lco_stat_fpath,"rb")) == NULL ) err(errno,"grouping_genomes():cannot find %s under %s ",sketch_stat,set_opt->insketchpath);
	dim_sketch_stat_t  lco_stat_readin;
  fread( &lco_stat_readin, sizeof(lco_stat_readin),1,tmpfh );
  fclose(tmpfh);

	if(lco_stat_readin.infile_num != subset->gn)
    err(errno,"grouping_genomes():%s's genome number %d not matches %s's genome number %d",lco_stat_fpath,lco_stat_readin.infile_num,set_opt->subsetf,subset->gn);
  free(lco_stat_fpath);
  mkdir_p(set_opt->outdir);
//read index
  char tmppath[PATHLEN]; 
	sprintf(tmppath,"%s/%s",set_opt->insketchpath,idx_sketch_suffix);
  if((tmpfh = fopen(tmppath,"rb")) == NULL)  err(errno,"lgrouping_genomes():%s",tmppath);
	uint64_t *tmp_idx  = malloc(sizeof(uint64_t) * (lco_stat_readin.infile_num + 1) ) ;
	fread(tmp_idx,sizeof(uint64_t), lco_stat_readin.infile_num + 1 ,tmpfh);
	fclose(tmpfh);
//read comblco
	sprintf(tmppath,"%s/%s",set_opt->insketchpath,combined_sketch_suffix);
	if((tmpfh = fopen(tmppath,"rb")) == NULL)  err(errno,"lgrouping_genomes():%s",tmppath);
	uint64_t *mem_comblco = malloc(sizeof(uint64_t) * tmp_idx[lco_stat_readin.infile_num]);
	fread(mem_comblco,sizeof(uint64_t),tmp_idx[lco_stat_readin.infile_num],tmpfh );
	fclose(tmpfh);

// out index
	int outfn = 0;
	uint64_t *out_idx  = calloc( (subset->taxn+1), sizeof(uint64_t));
//outf 
	sprintf(tmppath,"%s/%s",set_opt->outdir,combined_sketch_suffix);	
  if((tmpfh = fopen(tmppath,"wb")) == NULL)  err(errno,"lgrouping_genomes():%s",tmppath);	

	for(int t = 0; t < subset->taxn; t++){
		if(subset->tax[t].taxid == 0) continue;// ignore taxid 0		
		outfn++;
		out_idx[outfn] = out_idx[outfn-1];

		khash_t(kmer_set) *h = kh_init(kmer_set); int ret;		
		for(int n = 1; n <= subset->tax[t].gids[0];n++){
			int gid = subset->tax[t].gids[n] ;		
			for(uint64_t i = tmp_idx[gid]; i < tmp_idx[gid+1]; i++)
      	kh_put(kmer_set, h, mem_comblco[i],&ret);		
		}
		//write into comblco for tax t
		for (khint_t k = kh_begin(h); k != kh_end(h); ++k) {
      if (kh_exist(h, k)) {
				fwrite(&kh_key(h, k),sizeof(uint64_t),1,tmpfh);
				out_idx[outfn]++;	
			}
    }	
		kh_destroy(kmer_set,h);		
	}
	fclose(tmpfh);
	free(mem_comblco);
	free(tmp_idx);
	//write index
	sprintf(tmppath,"%s/%s",set_opt->outdir,idx_sketch_suffix);
	if((tmpfh = fopen(tmppath,"wb")) == NULL)  err(errno,"lgrouping_genomes():%s",tmppath);
	fwrite(out_idx,sizeof(uint64_t), outfn+1, tmpfh);
	free(out_idx); fclose(tmpfh);
	//write stat file 
  sprintf(tmppath,"%s/%s",set_opt->outdir,sketch_stat);
  if((tmpfh = fopen(tmppath,"wb")) == NULL)  err(errno,"lgrouping_genomes():%s",tmppath);
	lco_stat_readin.infile_num = outfn; lco_stat_readin.koc = 0;
	fwrite(&lco_stat_readin,sizeof(lco_stat_readin),1,tmpfh);

	char (*tmpfname)[PATHLEN] = malloc(outfn * PATHLEN); int idx = 0;
  for(int t=0;t<subset->taxn;t++){
    if(subset->tax[t].taxid != 0) {
      if(subset->tax[t].taxname != NULL)
        sprintf(tmpfname[idx],"%d_%s",subset->tax[t].taxid,subset->tax[t].taxname);
      else sprintf(tmpfname[idx],"%d",subset->tax[t].taxid);
      idx++;
    }
    free(subset->tax[t].gids);
    free(subset->tax[t].taxname);
  }
  free(subset->tax);
  free(subset);

  fwrite(tmpfname,PATHLEN,outfn,tmpfh);
	free(tmpfname);
  fclose(tmpfh);

	return outfn;
}

































