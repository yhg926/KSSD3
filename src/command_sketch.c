#include "command_sketch.h"
#include <time.h>

//shared global vars
const char sketch_stat[] = "lcofiles.stat";
const char sketch_suffix[] = "lco"; // long co, 64bits
const char combined_sketch_suffix[] = "comblco";
const char idx_sketch_suffix[] = "comblco.index";
const char combined_ab_suffix[] = "comblco.a";
//shared vars in this scope
static char TL;
static unsigned int FILTER;
static uint64_t ctxmask,tupmask;
//tmp container in this scope
static size_t file_size;
static dim_sketch_stat_t comblco_stat_one, comblco_stat_it;
static char tmp_fname[PATHLEN+20];
static struct stat tmpstat;
//functions
void compute_sketch(sketch_opt_t * sketch_opt_val, infile_tab_t* infile_stat){
	//initilization
	TL = 2*(sketch_opt_val->hclen + sketch_opt_val->holen) + sketch_opt_val->iolen;
	uint64_t tmp_var =  UINT64_MAX >> (64 - 2*sketch_opt_val->hclen) ;// (1LLU << (2*sketch_opt_val->hclen)) - 1 ;
	ctxmask = (tmp_var << (2*(sketch_opt_val->iolen + sketch_opt_val->hclen + sketch_opt_val->holen)) ) | (tmp_var << (2*(sketch_opt_val->holen) ));
 	tupmask = UINT64_MAX >> (64 - 2*TL) ;
	FILTER = UINT32_MAX >> sketch_opt_val->drfold  ; //2^(32-12)
	if ( TL > 32 || FILTER < 256) err(EINVAL,"compute_sketch(): TL (%d) or FILTER (%u) is out of range (TL <=32 and FILTER: 256..0xffffffff)",TL, FILTER);
	unsigned int hash_id =  GET_SKETCHING_ID(sketch_opt_val->hclen, sketch_opt_val->holen, sketch_opt_val->iolen, sketch_opt_val->drfold , FILTER); 
	printf("Sketching method hashid = %u\tTL=%u\n", hash_id,TL);
	uint64_t * tmp_ct_list = calloc (infile_stat->infile_num + 1, sizeof(uint64_t));	
	mkdir_p(sketch_opt_val->outdir);
#pragma omp parallel for num_threads(sketch_opt_val->p) schedule(guided)	
	for(int i = 0; i< infile_stat->infile_num; i++ ){
		// sketching all genomes individually
	 tmp_ct_list[i+1] = reads2sketch64( (infile_stat->organized_infile_tab)[i].fpath, format_string("%s/%d.%s",sketch_opt_val->outdir,i,sketch_suffix), sketch_opt_val->abundance, sketch_opt_val->kmerocrs);
		printf("\r%dth/%d genome sketching %s completed!\tsketch size=%lu",i+1,infile_stat->infile_num,(infile_stat->organized_infile_tab)[i].fpath,tmp_ct_list[i+1]);
	}
	printf("\n");
//combine *.lco to comblco
	FILE *comb_sketch_fp,*comb_ab_fp;
	if( ( comb_sketch_fp = fopen(format_string("%s/%s",sketch_opt_val->outdir,combined_sketch_suffix),"wb")) == NULL )  
		err(errno,"%s() open file error: %s/%s",__func__,sketch_opt_val->outdir,combined_sketch_suffix);	
 	if(sketch_opt_val->abundance){
		if( ( comb_ab_fp = fopen(format_string("%s/%s.a",sketch_opt_val->outdir,combined_sketch_suffix),"wb")) == NULL ) 
			err(errno,"%s() open file error: %s/%s.a",__func__,sketch_opt_val->outdir,combined_sketch_suffix);
	}
	
	for(int i = 0; i< infile_stat->infile_num; i++ ){	
		sprintf(tmp_fname,"%s/%d.%s",sketch_opt_val->outdir,i,sketch_suffix);
		uint64_t *mem_lco = read_from_file(tmp_fname, &file_size); 		
		fwrite(mem_lco,file_size,1,comb_sketch_fp);
		tmp_ct_list[i+1] += tmp_ct_list[i];	
		remove(tmp_fname);free(mem_lco);	
		//abundance 
		if(!sketch_opt_val->abundance) continue;
		sprintf(tmp_fname,"%s/%d.%s.a",sketch_opt_val->outdir,i,sketch_suffix);
		uint32_t *mem_ab = read_from_file(tmp_fname, &file_size);
		fwrite(mem_ab,file_size,1,comb_ab_fp);
		remove(tmp_fname);free(mem_ab);		
	}
	fclose(comb_sketch_fp);	if( sketch_opt_val->abundance) fclose(comb_ab_fp);	
 //write comblco.index 
	write_to_file(test_create_fullpath(sketch_opt_val->outdir,idx_sketch_suffix),tmp_ct_list,sizeof(tmp_ct_list[0]) * (infile_stat->infile_num + 1) );
  //write stat file
	dim_sketch_stat_t dim_sketch_stat = {
    .hash_id =  hash_id,
    .koc = sketch_opt_val->abundance,
    .klen = TL,
    .hclen = sketch_opt_val->hclen,
    .holen = sketch_opt_val->holen,
    .drfold = sketch_opt_val->drfold, //2^12 = 4096
    .infile_num = infile_stat->infile_num
  };
  char (*tmpname)[PATHLEN] = malloc(infile_stat->infile_num * PATHLEN);
  for(int i = 0; i< infile_stat->infile_num; i++ ) memcpy(tmpname[i],(infile_stat->organized_infile_tab)[i].fpath,PATHLEN);	
  concat_and_write_to_file(test_create_fullpath(sketch_opt_val->outdir,sketch_stat),&dim_sketch_stat,sizeof(dim_sketch_stat),tmpname,infile_stat->infile_num * PATHLEN );

	free_all( tmpname,tmp_ct_list,NULL);
};


KSEQ_INIT(gzFile, gzread) 
KHASH_MAP_INIT_INT64(kmer_hash, int)

int reads2sketch64 (char* seqfname, char * outfname, bool abundance, int n ) {
	
	gzFile infile = gzopen(seqfname, "r");	if( !infile ) err(errno,"reads2sketch64(): Cannot open file %s", seqfname);	
 	kseq_t *seq = kseq_init(infile);
	khash_t(kmer_hash) *h = kh_init(kmer_hash);
	uint64_t tuple,crvstuple,unituple,basenum,unictx;	
	uint32_t len_mv = 2*TL - 2;	 uint32_t sketch_size = 0;
 
	while (kseq_read(seq) >= 0) {
		const char *s = seq->seq.s;	 	
		if(seq->seq.l < TL) continue;	
		int base = 0;

		for(int pos = 0; pos < seq->seq.l ; pos++){ //for(int pos = TL; pos < seq->seq.l ; pos++)
			if (Basemap[(unsigned short)s[pos]] == DEFAULT){ base = 0;continue;}
			basenum = Basemap[(unsigned short)s[pos]];
      		tuple = ( ( tuple<< 2 ) | basenum )  ;
      		crvstuple = (( crvstuple >> 2 ) | ((basenum^3LLU) << len_mv )) ;
			if(++base < TL) continue;
			// if base >=TL, namely, contiue ACGT TL-mer
			unituple = (tuple & ctxmask) < (crvstuple & ctxmask) ? tuple : crvstuple;
      		unictx = unituple & ctxmask;
			if (SKETCH_HASH( unictx ) > FILTER) continue;						
			int ret; khint_t key = kh_put(kmer_hash, h, unituple & tupmask, &ret);

			if (ret) { 
				kh_value(h, key) = 1;
			 	sketch_size++;				
			} else {kh_value(h, key)++;}	
		}//for line		
	};//while
	kseq_destroy(seq);	gzclose(infile);

	//write sketch and abundance
	uint64_t *mem_lco = malloc(sketch_size * sizeof(uint64_t));
	uint32_t *mem_ab; 	if (abundance) mem_ab = malloc(sketch_size * sizeof(uint32_t));
	uint32_t kmer_ct = 0;
	for (khint_t k = kh_begin(h); k != kh_end(h); ++k) {
  	if (kh_exist(h, k) && kh_value(h, k) >= n ) {
			mem_lco[kmer_ct] = kh_key(h, k);
			if(abundance) mem_ab[kmer_ct]	= kh_value(h, k);
			kmer_ct++;	
		}
  }
	
	write_to_file(outfname,mem_lco,sketch_size * sizeof(uint64_t));
	kh_destroy(kmer_hash, h); free(mem_lco); 
	if(abundance){
		write_to_file(format_string("%s.a",outfname),mem_ab,sketch_size * sizeof(uint32_t));
		free(mem_ab);
	} 	
	return kmer_ct;// sketch_size;
}



int merge_comblco (sketch_opt_t * sketch_opt_val){

	void *mem_stat = read_from_file( test_get_fullpath(sketch_opt_val->remaining_args[0] ,sketch_stat), &file_size);
	memcpy(&comblco_stat_one,mem_stat,sizeof(comblco_stat_one));  
	comblco_stat_one.infile_num = 0; 
	char (*tmpname)[PATHLEN] = malloc(PATHLEN);
	FILE * merge_out_fp =  fopen( test_create_fullpath(sketch_opt_val->outdir,combined_sketch_suffix),"wb");
	if ( merge_out_fp == NULL ) err(errno,"%s():%s/%s",__func__,sketch_opt_val->outdir,combined_sketch_suffix);
	
	uint64_t *index_arry = malloc(sizeof(uint64_t)*100) ; index_arry[0] = 0;
	for(int i=0; i<sketch_opt_val->num_remaining_args; i++){
// read stat
		void *mem_stat_it = read_from_file( test_get_fullpath(sketch_opt_val->remaining_args[i] ,sketch_stat), &file_size);
		memcpy(&comblco_stat_it,mem_stat_it,sizeof(comblco_stat_it));
		if (comblco_stat_it.hash_id != comblco_stat_one.hash_id) err(EINVAL,"%uth %s hashid: %u != %u ",i,sketch_opt_val->remaining_args[i],comblco_stat_it.hash_id,comblco_stat_one.hash_id );
		if(comblco_stat_it.koc == 0 ) comblco_stat_one.koc = 0;
    //collect filenames
		tmpname = realloc(tmpname, PATHLEN * (comblco_stat_one.infile_num + comblco_stat_it.infile_num) );
		memcpy(tmpname+comblco_stat_one.infile_num,mem_stat_it+sizeof(comblco_stat_it),PATHLEN*comblco_stat_it.infile_num);
		//set index 
		uint64_t *mem_index_it = read_from_file( test_get_fullpath(sketch_opt_val->remaining_args[i],idx_sketch_suffix), &file_size);
		index_arry = (uint64_t*)realloc(index_arry,sizeof(uint64_t)*( comblco_stat_one.infile_num + comblco_stat_it.infile_num + 1));
		for(int j=1; j< comblco_stat_it.infile_num + 1;j++ ) index_arry[comblco_stat_one.infile_num + j] = index_arry[comblco_stat_one.infile_num] + mem_index_it[j];		
		// add file num
    comblco_stat_one.infile_num += comblco_stat_it.infile_num;		
		//write combined_sketch_suffix			
		uint64_t *mem_comblco = read_from_file( test_get_fullpath(sketch_opt_val->remaining_args[i] ,combined_sketch_suffix), &file_size) ;
		fwrite(mem_comblco,file_size, 1, merge_out_fp);			
		free_all(mem_stat_it,mem_index_it,mem_comblco,NULL);		
	}
	fclose(merge_out_fp); //write combined_sketch_suffix complete
	
	if(comblco_stat_it.koc){
		merge_out_fp =  fopen( test_create_fullpath(sketch_opt_val->outdir,combined_ab_suffix),"wb");
		if( merge_out_fp == NULL ) err(errno,"%s():%s/%s", __func__,sketch_opt_val->outdir,combined_ab_suffix);
		for(int i=0; i<sketch_opt_val->num_remaining_args; i++){
			uint64_t *mem_ab = read_from_file( test_get_fullpath(sketch_opt_val->remaining_args[i] ,combined_ab_suffix), &file_size) ;
			fwrite(mem_ab,file_size, 1, merge_out_fp);
			free(mem_ab);
		}
		fclose(merge_out_fp);	
	}
//write index and stat file
	write_to_file(test_create_fullpath(sketch_opt_val->outdir,idx_sketch_suffix), index_arry ,sizeof(index_arry[0])*(comblco_stat_one.infile_num+1));
	concat_and_write_to_file(test_create_fullpath(sketch_opt_val->outdir,sketch_stat),&comblco_stat_one,sizeof(comblco_stat_one),tmpname,PATHLEN * (comblco_stat_one.infile_num))	;

	return comblco_stat_one.infile_num;
}

































