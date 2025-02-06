#include "command_composite.h"
#include "global_basic.h"
#include "klib/khash.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <err.h>
#include <errno.h>
#include <math.h>
#include <libgen.h>
#include <dirent.h>

const char binVec_suffix[] = "abv";
const char abunMtx_suffix[] = "abm";
const char abunMtx_idx_suffix[] = "abmi";
const char abunMtx_name_suffix[] = "name";
const char binVec_dirname[] = "abundance_Vec";
const char y_l2n_suffix[] = "yl2n";

int read_abv (composite_opt_t *composite_opt){
	for(int i =0; i<composite_opt->num_remaining_args;i++){
    	const char *ext = strrchr(composite_opt->remaining_args[i],'.');
    	if(strcmp(ext + 1, binVec_suffix) != 0){
      		printf("%dth argument %s is not a .abv file, skipped\n",i,composite_opt->remaining_args[i]);
      		continue;
    	}
		size_t file_size; 
		binVec_t *tmp_abv = (binVec_t *)read_from_file( composite_opt->remaining_args[i], &file_size) ;
		int n = file_size/sizeof(binVec_t);
		printf("%d\t%s\n",n,composite_opt->remaining_args[i]);
		for(int l= 0 ; l < n; l++) printf("%d\t%f\n",tmp_abv[l].ref_idx,tmp_abv[l].pct);
		free(tmp_abv);
	}
	return 1;
}

float *tmp_measure;
typedef struct xny_pct { float x; float y;} xny_pct_t ; 

int abv_search	(composite_opt_t *composite_opt){		
		
	int abv_fn; size_t file_size;
	//read sample names, y_l2n,abunMtx idx and abunMtx
	char **tmpfname = read_lines_from_file( format_string("%s/%s.%s",composite_opt->refdir,binVec_dirname,abunMtx_name_suffix), &abv_fn);
	double *y_l2n = read_from_file(format_string("%s/%s.%s",composite_opt->refdir,binVec_dirname,y_l2n_suffix),&file_size);
	int *abunMtx_idx =  read_from_file(format_string("%s/%s.%s",composite_opt->refdir,binVec_dirname,abunMtx_idx_suffix),&file_size);
	binVec_t *abunMtx = read_from_file(format_string("%s/%s.%s",composite_opt->refdir,binVec_dirname,abunMtx_suffix),&file_size);
	
	xny_pct_t *xny_pct ;// for l1norm calculation
	int *abv_ids = malloc(sizeof(int) * abv_fn);// matched abv ids
#define DFLT (-2) // set default measure, value should not in the range of l1norm l2norm cosine   	
	tmp_measure = malloc(sizeof(float) * abv_fn);
	if(composite_opt->s == 1) xny_pct = malloc(sizeof(xny_pct_t) *abv_fn);

	for(int i =0; i<composite_opt->num_remaining_args;i++){
		const char *ext = strrchr(composite_opt->remaining_args[i],'.');
		if(strcmp(ext + 1, binVec_suffix) != 0){
			printf("%dth argument %s is not a .abv file, skipped\n",i,composite_opt->remaining_args[i]);
			continue;
		}
		int n_abv_match = 0; //num of matched abv
		char *abv_fpath;  // absolute path
		if( strchr(composite_opt->remaining_args[i],'/') != NULL ) abv_fpath = composite_opt->remaining_args[i];
		else abv_fpath = format_string("%s/%s/%s", composite_opt->refdir,binVec_dirname,composite_opt->remaining_args[i]);
		binVec_t *tmp_abv = (binVec_t *)read_from_file(abv_fpath, &file_size);
		int abv_len = file_size/sizeof(binVec_t);
	
		for(int i = 0; i<abv_fn ; i++) tmp_measure[i] = DFLT;
		if(composite_opt->s == 1) memset(xny_pct,0, sizeof(xny_pct_t) *abv_fn );
		
		float xl2n = 0; //query l2norm ||x|| 
		for(int d = 0; d< abv_len; d++){
			int ref_idx = tmp_abv[d].ref_idx;
			xl2n += tmp_abv[d].pct*tmp_abv[d].pct;
			for(int j = ref_idx==0?0:abunMtx_idx[ref_idx-1]; j< abunMtx_idx[ref_idx]; j++){

				if(tmp_measure[abunMtx[j].ref_idx]  == DFLT) {
					tmp_measure[abunMtx[j].ref_idx] = 0;
					abv_ids[n_abv_match] = abunMtx[j].ref_idx;
					n_abv_match++;
				}
				if(composite_opt->s == 1){
					tmp_measure[abunMtx[j].ref_idx] += (float)fabs(abunMtx[j].pct - tmp_abv[d].pct);
					xny_pct[abunMtx[j].ref_idx].x += tmp_abv[d].pct;
					xny_pct[abunMtx[j].ref_idx].y += abunMtx[j].pct;
				}
				else if (composite_opt->s == 2)
					tmp_measure[abunMtx[j].ref_idx] += (abunMtx[j].pct - tmp_abv[d].pct) * (abunMtx[j].pct - tmp_abv[d].pct);
				else //cosine
					tmp_measure[abunMtx[j].ref_idx] += abunMtx[j].pct * tmp_abv[d].pct;					
			}
		}//dimension
		
		if(composite_opt->s == 0) {
			for(int n = 0; n < n_abv_match; n++)
				tmp_measure[abv_ids[n]] = tmp_measure[abv_ids[n]] / (sqrt(xl2n)*y_l2n[abv_ids[n]]) ;
		}		
		printf("#Sample\t");
		//sort
		if(composite_opt->s == 1) {
			for(int n = 0; n < n_abv_match; n++) 
				tmp_measure[abv_ids[n]] += (2*100 - xny_pct[abv_ids[n]].x - xny_pct[abv_ids[n]].y); // 2*100 bcs abv has been scaled to 100

			qsort(abv_ids,n_abv_match,sizeof(int),comparator_measure);
			printf("L1norm\n");
			for(int n = 0; n < n_abv_match; n++) printf("%s\t%lf\n",tmpfname[abv_ids[n]],tmp_measure[abv_ids[n]]);
		}
		else if(composite_opt->s == 2){
			qsort(abv_ids,n_abv_match,sizeof(int),comparator_measure);
			printf("L2norm\n");
			for(int n = 0; n < n_abv_match; n++) printf("%s\t%lf\n",tmpfname[abv_ids[n]],sqrt(tmp_measure[abv_ids[n]]));
		}
		else{
			qsort(abv_ids,n_abv_match,sizeof(int),comparator_measure);
			printf("CosineXY\n");
			for(int n = n_abv_match -1 ; n >= 0; n--) printf("%s\t%lf\n",tmpfname[abv_ids[n]],tmp_measure[abv_ids[n]]);	
		}

		free(tmp_abv);
	}//for remaining_args[i]

	free_all(tmp_measure,tmpfname,abunMtx,abunMtx_idx,y_l2n,NULL);
	if(composite_opt->s == 1) free(xny_pct);

	return 1;
}


int index_abv (composite_opt_t *composite_opt){

	char* abv_dpath = format_string("%s/%s", composite_opt->refdir,binVec_dirname);
	DIR *abv_dir = opendir(abv_dpath); if(!abv_dir) err(errno, "index_abv(): %s does not exists! \n",abv_dpath);
	unify_sketch_t *result = generic_sketch_parse(composite_opt->refdir);
	int infile_num = result->infile_num; 
	free_unify_sketch(result); 

	int *abunMtx_count = calloc(infile_num,sizeof(int));
	binVec_t **abunMtx_dynamic = malloc(sizeof(binVec_t *) * infile_num);
	for(int i =0; i< infile_num;i++) abunMtx_dynamic[i] = malloc(sizeof(binVec_t));	

	int abv_fcnt = 0; size_t file_size;
	FILE *namefh,*abunMtx_fh,*y_l2n_fh;
	struct dirent *dirent; binVec_t binVec_tmp;	
	char* abv_fpath = malloc(PATHLEN);
	sprintf(abv_fpath,"%s/%s.%s",composite_opt->refdir,binVec_dirname,abunMtx_name_suffix);
	if(( namefh = fopen(abv_fpath,"w") ) == NULL) err(errno, "index_abv():%s",abv_fpath);
	sprintf(abv_fpath,"%s/%s.%s",composite_opt->refdir,binVec_dirname,y_l2n_suffix);
	if(( y_l2n_fh = fopen(abv_fpath,"wb") ) == NULL) err(errno, "index_abv():%s",abv_fpath);

	while (( dirent = readdir(abv_dir)) != NULL){
		if(strcmp(strrchr(dirent->d_name,'.') + 1, binVec_suffix) != 0) continue;
		binVec_t *binVec_tmp = (binVec_t *)read_from_file( format_string("%s/%s",abv_dpath,dirent->d_name), &file_size); 
		double y_l2n = 0; int abv_len = file_size/sizeof(binVec_t) ;
		for(int i =0; i<abv_len ;i++){
			y_l2n+= binVec_tmp[i].pct*binVec_tmp[i].pct;
			abunMtx_dynamic[binVec_tmp[i].ref_idx] = realloc((binVec_t *)abunMtx_dynamic[binVec_tmp[i].ref_idx], sizeof(binVec_t) * (abunMtx_count[binVec_tmp[i].ref_idx] + 1) );
       		abunMtx_dynamic[binVec_tmp[i].ref_idx][abunMtx_count[binVec_tmp[i].ref_idx]].ref_idx = abv_fcnt ;
			abunMtx_dynamic[binVec_tmp[i].ref_idx][abunMtx_count[binVec_tmp[i].ref_idx]].pct = binVec_tmp[i].pct;
			abunMtx_count[binVec_tmp[i].ref_idx]++;
		}
		free(binVec_tmp);						

		sprintf(abv_fpath,"%s\n",dirent->d_name);
		fwrite(abv_fpath,strlen(abv_fpath),1,namefh);			
		y_l2n = sqrt(y_l2n);
		fwrite(&y_l2n,sizeof(double),1,y_l2n_fh);
		abv_fcnt++;							
	};
	closedir(abv_dir); fclose(namefh); fclose(y_l2n_fh);

 	sprintf(abv_fpath,"%s/%s.%s",composite_opt->refdir,binVec_dirname,abunMtx_suffix);
	if(( abunMtx_fh = fopen(abv_fpath,"wb") ) == NULL) err(errno, "index_abv():%s",abv_fpath);
	for(int i = 0 ;i < infile_num;i++){
		if(abunMtx_count[i])	
			fwrite( abunMtx_dynamic[i], sizeof(binVec_t), abunMtx_count[i],abunMtx_fh);
	}
	fclose(abunMtx_fh);

	for(int i = 1 ;i < infile_num;i++) abunMtx_count[i] += abunMtx_count[i-1];
	write_to_file(format_string("%s/%s.%s",composite_opt->refdir,binVec_dirname,abunMtx_idx_suffix),abunMtx_count,sizeof(int)*infile_num );

	free_all(abv_dpath,abv_fpath,abunMtx_count,abunMtx_dynamic,NULL);

	return 1;
}

// this version get_species_abundance () assume few qry input genome,
// it waste time in reading qry in for qry.infile_num loop 
KHASH_MAP_INIT_INT64(kmer_hash, int)

#ifndef MIN_KM_S
#define MIN_KM_S 6 // minimal kmer share for a ref genome
#endif
#define ST_PCTL (0.98) //start percentile
#define ED_PCTL (0.99) //end percentile
Vector *ref_abund;  //make it global, seen by qsort comparetor
int get_species_abundance (composite_opt_t * composite_opt) { //by uniq kmer in a species

    unify_sketch_t* ref_result = generic_sketch_parse(composite_opt->refdir), *qry_result = generic_sketch_parse(composite_opt->qrydir);
    if(ref_result->stat_type != qry_result->stat_type) err(EXIT_FAILURE,"%s(): ref sketch type %u != qry %u",__func__,ref_result->stat_type,qry_result->stat_type);
    else if(ref_result->hash_id != qry_result->hash_id) err(EXIT_FAILURE,"%s(): ref hash_id %u != qry %u",__func__,ref_result->hash_id,qry_result->hash_id);	
	if(qry_result->abundance == NULL) err(errno, "%s():query has no abundance", __func__);	
	//initialize sort ref, will sort gid by shared k-mer for each qry 
    int *sort_ref = malloc(ref_result->infile_num* sizeof(int));
    for(int i = 0; i< ref_result->infile_num; i++) sort_ref[i] = i;
	ref_abund = malloc(ref_result->infile_num* sizeof(Vector));
	binVec_t* binVec = malloc(ref_result->infile_num* sizeof(binVec_t));

	for(int qn = 0; qn < qry_result->infile_num; qn++) { // for each query genome, not suite for qry num is large
		//initialize hashtalbe for each qry genome
		khash_t(kmer_hash) *km2abund = kh_init(kmer_hash); int ret;
		for (uint64_t idx = qry_result->sketch_index[qn] ; idx < qry_result->sketch_index[qn+1]; idx++){
            khint_t key = kh_put(kmer_hash, km2abund, qry_result->comb_sketch[idx], &ret);
			kh_value(km2abund, key) = qry_result->abundance[idx] ;
		}
#pragma omp parallel for num_threads(composite_opt->p) schedule(guided)
		for(int rn = 0; rn < ref_result->infile_num; rn++) { // for each ref genome	
			vector_init(&ref_abund[rn], sizeof(int));//initialize vector
			for( uint64_t ri = ref_result->sketch_index[rn]; ri < ref_result->sketch_index[rn+1]; ri++){
				khint_t key = kh_get(kmer_hash, km2abund, ref_result->comb_sketch[ri] ) ; 
				if(key != kh_end(km2abund)) vector_push( &ref_abund[rn], &kh_value(km2abund, key));						
			}		
		}
		kh_destroy(kmer_hash, km2abund);		
		// sort ref_abund
		qsort(sort_ref, ref_result->infile_num, sizeof(sort_ref[0]), comparator_idx);		
//if binary vector
//		binVec_t* binVec = malloc(ref_result->infile_num* sizeof(binVec_t));
		int num_pass = 0;		float binVecsum = 0;

		for ( int i = 0; i< ref_result->infile_num; i++ ){
			int kmer_num = ref_abund[sort_ref[i]].size; // overlapped kmer_num 
			if (kmer_num < MIN_KM_S) break; 						
			qsort(ref_abund[sort_ref[i]].data, kmer_num, sizeof(int),comparator);
//average
			int sum = 0; 
			for(int n = 0; n < kmer_num; n++) sum += *(int *)vector_get(&ref_abund[sort_ref[i]], n) ;// (int *)(ref_abund[sort_ref[i]].data) [n];
//median:real median is (kmer_num + 1)/2, but we ignore the last(most abundant) k-mer so use			
			int median_idx = kmer_num /2;
			int pct09_idx = kmer_num  * ST_PCTL ;
			int median = *(int *)vector_get(&ref_abund[sort_ref[i]], median_idx);
			int max = *(int *)vector_get(&ref_abund[sort_ref[i]], kmer_num -1 );  
//percentile range average
			int lastsum = 0; int lastn = kmer_num*ED_PCTL - pct09_idx + 1;
      		for(int n = pct09_idx ; n < kmer_num*ED_PCTL; n++) lastsum += *(int *)vector_get(&ref_abund[sort_ref[i]], n); //ref_abund[sort_ref[i]][n];
		
			vector_free(&ref_abund[sort_ref[i]]);				
		
			if(!composite_opt->b)
				printf("%s\t%s\t%d\t%f\t%f\t%d\t%d\n",qry_result->gname[qn],ref_result->gname[sort_ref[i]], kmer_num, (float)sum/kmer_num,(float)lastsum/lastn, median, max );
			else if ( median > 1 && kmer_num > MIN_KM_S + 1){ //  threshold for abv
				binVec[num_pass].ref_idx = sort_ref[i];
				binVec[num_pass].pct = (float)lastsum/lastn;
				binVecsum += binVec[num_pass].pct;
				num_pass++;
			}			
		}
//output binary vector
		if(composite_opt->b){
			for(int i = 0; i< num_pass;i++) binVec[i].pct = ( binVec[i].pct - 1)*100/(binVecsum - num_pass) ;		
			char tmpfname[PATHLEN]; 		
		 	if(strlen(composite_opt->outdir) < 3) sprintf(tmpfname,"%s/%s",composite_opt->refdir,binVec_dirname);// outdir not set, value:\"./\"
            else strcpy(tmpfname,composite_opt->outdir);// if -o set, abv will be in outdir, not binVec_dirname
            mkdir_p(tmpfname);
	        write_to_file(format_string("%s/%s.%s",tmpfname,basename(qry_result->gname[qn]),binVec_suffix),binVec,sizeof(binVec_t)*num_pass );
		}	
	}// for qry
  	free_all(sort_ref,ref_abund,binVec,NULL) ;
	free_unify_sketch(ref_result); free_unify_sketch(qry_result);
	return 1;
}

inline int comparator_idx (const void *a, const void *b){
	return ( ref_abund[*(int*)b].size - ref_abund[*(int*)a].size );		
}

inline int comparator (const void *a, const void *b){
  return ( *(int*)a - *(int*)b );
}

inline int comparator_measure (const void *a, const void *b){
	float rtv = tmp_measure[*(int*)a] - tmp_measure[*(int*)b] ; 
    if(rtv > 0) return 1;
	else if(rtv < 0) return -1;
	else return 0;
}

