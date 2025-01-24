#include "command_sketch.h"
#include <time.h>

//shared global vars
static char TL;
static unsigned int FILTER;
static llong ctxmask,tupmask;

const char sketch_stat[] = "lcofiles.stat"; 
const char sketch_suffix[] = "lco"; // long co, 64bits
const char combined_sketch_suffix[] = "comblco";
const char idx_sketch_suffix[] = "comblco.index";
const char combined_ab_suffix[] = "comblco.a";

void compute_sketch(sketch_opt_t * sketch_opt_val, infile_tab_t* infile_stat){
	//initilization
	TL = 2*(sketch_opt_val->hclen + sketch_opt_val->holen) + sketch_opt_val->iolen;
	FILTER = UINT32_MAX >> sketch_opt_val->drfold  ; //2^(32-12)
	if ( TL > 32 || FILTER < 256) err(EINVAL,"compute_sketch(): TL (%d) or FILTER (%u) is out of range (TL <=32 and FILTER: 256..0xffffffff)",TL, FILTER);
	unsigned int hash_id =  GET_SKETCHING_ID(sketch_opt_val->hclen, sketch_opt_val->holen, sketch_opt_val->iolen, sketch_opt_val->drfold , FILTER); 
	printf("Sketching method hashid = %u\tTL=%u\n", hash_id,TL);
	dim_sketch_stat_t dim_sketch_stat = {
    .hash_id =  hash_id,
		.koc = sketch_opt_val->abundance,
    .klen = TL,
    .hclen = sketch_opt_val->hclen,
    .holen = sketch_opt_val->holen,
    .drfold = sketch_opt_val->drfold, //2^12 = 4096
		.infile_num = infile_stat->infile_num
  };
	if(mkdir_p(sketch_opt_val->outdir) != 0)  err(EINVAL,"Failed to create output directories:%s",sketch_opt_val->outdir);
	//write stat file
	char *tmp_outfname = malloc(PATHLEN+10);
  sprintf(tmp_outfname,"%s/%s",sketch_opt_val->outdir,sketch_stat);
	FILE *sketch_stat_fp;  
	if( ( sketch_stat_fp = fopen(tmp_outfname,"wb")) == NULL )  err(errno,"compute_sketch() open file error: %s",tmp_outfname);
  fwrite(&dim_sketch_stat,sizeof(dim_sketch_stat),1,sketch_stat_fp);
	for(int i = 0; i< infile_stat->infile_num; i++ )
		fwrite((infile_stat->organized_infile_tab)[i].fpath,PATHLEN,1, sketch_stat_fp);
	fclose(sketch_stat_fp);
	
	llong tmp_var = (1LLU << (2*sketch_opt_val->hclen)) - 1 ;	
	ctxmask = (tmp_var << (2*(sketch_opt_val->iolen + sketch_opt_val->hclen + sketch_opt_val->holen)) )
         + (tmp_var << (2*(sketch_opt_val->holen) ));

	tupmask = UINT64_MAX >> (64 - 2*TL) ;  

	// sketching individually
	uint64_t * tmp_ct_list = calloc (infile_stat->infile_num + 1, sizeof(uint64_t));

	uint32_t max_sketch_size = 0; 
#pragma omp parallel for num_threads(sketch_opt_val->p) reduction(max:max_sketch_size) schedule(guided)	
	for(int i = 0; i< infile_stat->infile_num; i++ ){
		char *tmp_fname = malloc(PATHLEN+10);	
		sprintf(tmp_fname,"%s/%d.%s",sketch_opt_val->outdir,i,sketch_suffix);

	 	tmp_ct_list[i+1] = reads2sketch64( (infile_stat->organized_infile_tab)[i].fpath, tmp_fname, sketch_opt_val->abundance);
		
		if(tmp_ct_list[i+1] > max_sketch_size)  max_sketch_size = tmp_ct_list[i+1];
		free(tmp_fname);
		printf("sketching %s completed!\tsketch size=%lu\n",(infile_stat->organized_infile_tab)[i].fpath,tmp_ct_list[i+1]);
	}


//combine *.lco to comblco
	void *tmp_mem = malloc(max_sketch_size*sizeof(llong));
	struct stat tmpstat;
	sprintf(tmp_outfname,"%s/%s",sketch_opt_val->outdir,combined_sketch_suffix);
	FILE *comb_sketch_fp,*comb_ab_fp;
	if( ( comb_sketch_fp = fopen(tmp_outfname,"wb")) == NULL )  err(errno,"compute_sketch() open file error: %s",tmp_outfname);	
 	if(sketch_opt_val->abundance){
		sprintf(tmp_outfname,"%s/%s.a",sketch_opt_val->outdir,combined_sketch_suffix);
		if( ( comb_ab_fp = fopen(tmp_outfname,"wb")) == NULL )  
			err(errno,"compute_sketch() open file error: %s",tmp_outfname);
	}
	char *tmp_infname = malloc(PATHLEN+10);	
	for(int i = 0; i< infile_stat->infile_num; i++ ){	
		sprintf(tmp_infname,"%s/%d.%s",sketch_opt_val->outdir,i,sketch_suffix);
		stat(tmp_infname,&tmpstat);
		FILE *sketch_fp;
		if( ( sketch_fp = fopen(tmp_infname,"rb")) == NULL )  err(errno,"compute_sketch() open file error: %s",tmp_infname);
	
		fread(tmp_mem,tmpstat.st_size,1,sketch_fp);
		fwrite(tmp_mem,tmpstat.st_size,1,comb_sketch_fp);			 		
		tmp_ct_list[i+1] += tmp_ct_list[i];
		fclose(sketch_fp);
		remove(tmp_infname);	

		if(!sketch_opt_val->abundance) continue;
		sprintf(tmp_infname,"%s/%d.%s.a",sketch_opt_val->outdir,i,sketch_suffix);	
		stat(tmp_infname,&tmpstat);
		if( ( sketch_fp = fopen(tmp_infname,"rb")) == NULL )  err(errno,"compute_sketch() open file error: %s",tmp_infname);
		fread(tmp_mem,tmpstat.st_size,1,sketch_fp);
	  fwrite(tmp_mem,tmpstat.st_size,1,comb_ab_fp);
		fclose(sketch_fp);
		remove(tmp_infname);		
	}
	free(tmp_infname);
	free(tmp_mem);
	fclose(comb_sketch_fp);
	if( sketch_opt_val->abundance) fclose(comb_ab_fp);	

 //write comblco.index
	sprintf(tmp_outfname,"%s/%s",sketch_opt_val->outdir,idx_sketch_suffix);
	FILE *idx_sketch_fp;
	if( ( idx_sketch_fp = fopen(tmp_outfname,"wb")) == NULL )  err(errno,"compute_sketch() open file error: %s",tmp_outfname);		
	fwrite(tmp_ct_list,sizeof(tmp_ct_list[0]),infile_stat->infile_num + 1,idx_sketch_fp);
	fclose(idx_sketch_fp);
	free(tmp_ct_list);	
	free(tmp_outfname);

};


KSEQ_INIT(gzFile, gzread) 
KHASH_MAP_INIT_INT64(kmer_hash, int)

int reads2sketch64 (char* seqfname, char * outfname, bool abundance ) {
	
	gzFile infile = gzopen(seqfname, "r");
	if( !infile ) err(errno,"reads2sketch64(): Cannot open file %s", seqfname);
	
 	kseq_t *seq = kseq_init(infile);
	khash_t(kmer_hash) *h = kh_init(kmer_hash);

	llong tuple,crvstuple,unituple,basenum,unictx;	
	int len_mv = 2*TL - 2;	

	while (kseq_read(seq) >= 0) {
		const char *s = seq->seq.s;	 	
		if(seq->seq.l < TL) continue;	

		for(int pos = 0; pos < TL; pos++){
      basenum = (llong)Basemap[(unsigned short)s[pos]];
			tuple = ( ( tuple<< 2 ) | basenum )  ;
      crvstuple = (( crvstuple >> 2 ) | ((basenum^3LLU) << len_mv )) ;
		}		
		for(int pos = TL; pos < seq->seq.l ; pos++){
			basenum =	(llong)Basemap[(unsigned short)s[pos]];
      tuple = ( ( tuple<< 2 ) | basenum )  ;
      crvstuple = (( crvstuple >> 2 ) | ((basenum^3LLU) << len_mv )) ;

			unituple = (tuple & ctxmask) < (crvstuple & ctxmask) ? tuple : crvstuple;
      unictx = unituple & ctxmask;

			if (SKETCH_HASH( unictx ) > FILTER) continue;
				
			int ret;
			khint_t key = kh_put(kmer_hash, h, unituple & tupmask, &ret);
			if (ret) { kh_value(h, key) = 1;}
      else {kh_value(h, key)++;}	
		}//for line		
	};//while
	kseq_destroy(seq);
	gzclose(infile);
	//write sketch

	FILE *sketch_fp;
	if( (sketch_fp =  fopen(outfname,"wb")) == NULL )  err(errno,"reads2sketch64() open file error: %s",outfname);
	uint32_t sketch_size = 0; 

	for (khint_t k = kh_begin(h); k != kh_end(h); ++k) {
  	if (kh_exist(h, k)) {
    	fwrite(&kh_key(h, k),sizeof(llong),1,sketch_fp );
			sketch_size++;
		}
  }
	
	fclose(sketch_fp);
	//write abundnace

	if(abundance) {
		FILE *ab_fp; char ab_fname[PATHLEN];
		sprintf(ab_fname,"%s.a",outfname);
  	if( (ab_fp =  fopen(ab_fname,"wb")) == NULL )  err(errno,"reads2sketch64() open file error: %s",ab_fname);
  	for (khint_t k = kh_begin(h); k != kh_end(h); ++k) {
    	if (kh_exist(h, k))  fwrite(&kh_value(h, k), sizeof(uint32_t),1, ab_fp );
  	}	
  	fclose(ab_fp);
	}

	kh_destroy(kmer_hash, h);

	return sketch_size;
}


int merge_comblco (sketch_opt_t * sketch_opt_val){
	if(mkdir_p(sketch_opt_val->outdir) != 0)  err(EINVAL,"Failed to create output directories:%s",sketch_opt_val->outdir);	
	dim_sketch_stat_t comblco_stat_one, comblco_stat_it;
	FILE *merge_out_fp, *tmp_fp;
	char *tmp_fpath = test_get_fullpath(sketch_opt_val->remaining_args[0], sketch_stat);	
	if( ( tmp_fp = fopen(tmp_fpath,"rb")) == NULL ) err(errno,"merge_comblco():%s", tmp_fpath);
	free(tmp_fpath);

	fread(&comblco_stat_one, sizeof(dim_sketch_stat_t), 1, tmp_fp);
	comblco_stat_one.infile_num = 0;
	fclose(tmp_fp);

	char (*tmpname)[PATHLEN] = malloc(PATHLEN);
	void *tmp_mem =  malloc(100);
	uint64_t  tmp_arry[4];
	uint64_t *index_arry = malloc(sizeof(uint64_t)*100) ; index_arry[0] = 0;

	struct stat file_stat;
	
  char merge_out_fn[PATHLEN+20]; 
 	sprintf( merge_out_fn,"%s/%s",sketch_opt_val->outdir,combined_sketch_suffix);

  if( ( merge_out_fp = fopen(merge_out_fn,"wb")) == NULL ) err(errno,"merge_comblco():%s", merge_out_fn);

	for(int i=0; i<sketch_opt_val->num_remaining_args; i++){
		tmp_fpath = test_get_fullpath(sketch_opt_val->remaining_args[i], sketch_stat);
		if( ( tmp_fp = fopen(tmp_fpath,"rb")) == NULL ) err(errno,"merge_comblco():%s", tmp_fpath);
		free(tmp_fpath);
 
		fread(&comblco_stat_it, sizeof(dim_sketch_stat_t), 1, tmp_fp);
		if (comblco_stat_it.hash_id != comblco_stat_one.hash_id) err(EINVAL,"%uth %s hashid: %u != %u ",i,sketch_opt_val->remaining_args[i],comblco_stat_it.hash_id,comblco_stat_one.hash_id );
		if(comblco_stat_it.koc == 0 ) comblco_stat_one.koc = 0;

    //collect filenames
		tmpname = realloc(tmpname, PATHLEN * (comblco_stat_one.infile_num + comblco_stat_it.infile_num) ) ;
		fread(tmpname+comblco_stat_one.infile_num,PATHLEN,comblco_stat_it.infile_num,tmp_fp);
		fclose(tmp_fp);

		//set index 
		
		tmp_fpath = test_get_fullpath(sketch_opt_val->remaining_args[i], idx_sketch_suffix);
		if( ( tmp_fp = fopen(tmp_fpath,"rb")) == NULL ) err(errno,"merge_comblco(): %dth arg: %s has no %s", i,tmp_fpath,idx_sketch_suffix);
		free(tmp_fpath);

		tmp_mem = (uint64_t*)realloc(tmp_mem,sizeof(uint64_t)*(comblco_stat_it.infile_num + 1));
		fread(tmp_mem, sizeof(uint64_t), comblco_stat_it.infile_num + 1, tmp_fp);
		index_arry = (uint64_t*)realloc(index_arry,sizeof(uint64_t)*( comblco_stat_one.infile_num + comblco_stat_it.infile_num + 1));

		for(int j=1; j< comblco_stat_it.infile_num + 1;j++ )
			  	index_arry[comblco_stat_one.infile_num + j] = index_arry[comblco_stat_one.infile_num] + *((uint64_t*)tmp_mem+j);
		
		fclose(tmp_fp);

		// add file num
    comblco_stat_one.infile_num += comblco_stat_it.infile_num;		

		//write combined_sketch_suffix	
		tmp_fpath = test_get_fullpath(sketch_opt_val->remaining_args[i],combined_sketch_suffix);
		if( ( tmp_fp = fopen(tmp_fpath,"rb")) == NULL ) err(errno,"merge_comblco(): %dth arg: %s has no %s", i,tmp_fpath,combined_sketch_suffix);
		stat(tmp_fpath, &file_stat);
		tmp_mem = realloc(tmp_mem,file_stat.st_size);
		fread(tmp_mem, file_stat.st_size, 1, tmp_fp);
		fwrite(tmp_mem,file_stat.st_size, 1, merge_out_fp);
		fclose(tmp_fp);	free(tmp_fpath);					
	}
	fclose(merge_out_fp); //write combined_sketch_suffix complete
	
	if(comblco_stat_it.koc){

		sprintf( merge_out_fn,"%s/%s",sketch_opt_val->outdir,combined_ab_suffix);
		if( ( merge_out_fp = fopen(merge_out_fn,"wb")) == NULL ) err(errno,"merge_comblco():%s", merge_out_fn);

		for(int i=0; i<sketch_opt_val->num_remaining_args; i++){
			tmp_fpath = test_get_fullpath(sketch_opt_val->remaining_args[i],combined_ab_suffix);
			if( ( tmp_fp = fopen(tmp_fpath,"rb")) == NULL ) err(errno,"merge_comblco(): %dth arg: %s has no %s", i,tmp_fpath,combined_ab_suffix);
			stat(tmp_fpath, &file_stat);
			tmp_mem = realloc(tmp_mem,file_stat.st_size);
			fread(tmp_mem, file_stat.st_size, 1, tmp_fp);
			fwrite(tmp_mem,file_stat.st_size, 1, merge_out_fp);
			fclose(tmp_fp); free(tmp_fpath);			
		}
		fclose(merge_out_fp);	
	}
	free(tmp_mem);
//write index
	sprintf( merge_out_fn,"%s/%s",sketch_opt_val->outdir,idx_sketch_suffix);		
	if( ( merge_out_fp = fopen(merge_out_fn,"wb")) == NULL ) err(errno,"merge_comblco():%s", merge_out_fn);
	fwrite(index_arry,sizeof(uint64_t),comblco_stat_one.infile_num+1,merge_out_fp);
	fclose(merge_out_fp);
	free(index_arry);

//write stat file
	sprintf( merge_out_fn,"%s/%s",sketch_opt_val->outdir,sketch_stat);
	if( ( merge_out_fp = fopen(merge_out_fn,"wb")) == NULL ) err(errno,"merge_comblco():%s", merge_out_fn);
	fwrite(&comblco_stat_one,sizeof(dim_sketch_stat_t), 1,merge_out_fp);
	fwrite(tmpname,PATHLEN,comblco_stat_one.infile_num,merge_out_fp);	
	free(tmpname);
	fclose(merge_out_fp);
	return comblco_stat_one.infile_num;
}

































