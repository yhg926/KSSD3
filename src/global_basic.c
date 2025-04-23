/*used by both global and sumcommand*/ 
#include "global_basic.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <err.h>
#include <errno.h>
#include <ctype.h>  
#include <sys/stat.h>
#include <dirent.h>
#include <sys/sysinfo.h>
#include <math.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>


/*global basemap and mapbase*/
#if ALPHABET == 1 // objs compatible mode
const int Basemap[128] =
{ //make sure its ok illegal base letter has same behavial with null base
  [0 ... 127] = DEFAULT, ['z'] = DEFAULT, ['Z'] = DEFAULT ,
  ['a'] = 0, ['A'] = 0, ['c'] = 1, ['C'] = 1, ['g'] = 2, ['G'] = 2, ['t'] = 3, ['T'] = 3,
  ['w'] = 4, ['W'] = 4, ['s'] = 5, ['S'] = 5, ['m'] = 6, ['M'] = 6, ['k'] = 7, ['K'] = 7,
  ['r'] = 8, ['R'] = 8, ['y'] = 9, ['Y'] = 9, ['b'] = 10, ['B'] = 10, ['d'] = 11, ['D'] = 11,
  ['h'] = 12, ['H'] = 12, ['v'] = 13, ['V'] = 13, ['n'] = 14, ['N'] = 14,
};
const char Mapbase[] = {'A','C','G','T','W','S','M','K','R','Y','B','D','H','V','N','Z'};

#include <stdbool.h>
const bool Objdist[16][16] =
{
  [0 ... 15][0 ... 15] = 0,

  [0][1] = 1, [0][2] = 1, [0][3] = 1, [0][5] = 1, [0][7] = 1, [0][9] = 1, [0][10] = 1,
  [1][0] = 1, [2][0] = 1, [3][0] = 1, [5][0] = 1, [7][0] = 1, [9][0] = 1, [10][0] = 1,

  [1][2] = 1, [1][3] = 1, [1][4] = 1, [1][7] = 1, [1][8] = 1, [1][11] = 1,
  [2][1] = 1, [3][1] = 1, [4][1] = 1, [7][1] = 1, [8][1] = 1, [11][1] = 1,

  [2][3] = 1, [2][4] = 1, [2][6] = 1, [2][9] = 1, [2][12] = 1,
  [3][2] = 1, [4][2] = 1, [6][2] = 1, [9][2] = 1, [12][2] = 1,

  [3][5] = 1, [3][6] = 1, [3][8] = 1, [3][13] = 1,
  [5][3] = 1, [6][3] = 1, [8][3] = 1, [13][3] = 1,

  [4][5] = 1,   [6][7] = 1, [8][9] = 1,
  [5][4] = 1,   [7][6] = 1, [9][8] = 1,
};

#elif ALPHABET == 2 // AA sequences

const int Basemap[128] =
{
  [0 ... 127] = DEFAULT, ['a'] = 0, ['A'] = 0, ['c'] = 1, ['C'] = 1, ['d'] = 2, ['D'] = 2, ['e'] = 3, ['E'] = 3,
  ['f'] = 4, ['F'] = 4, ['g'] = 5, ['G'] = 5, ['h'] = 6, ['H'] = 6, ['i'] = 7, ['I'] = 7,
  ['k'] = 8, ['K'] = 8, ['l'] = 9, ['L'] = 9, ['m'] = 10, ['M'] = 10, ['n'] = 11, ['N'] = 11,
  ['p'] = 12, ['P'] = 12, ['q'] = 13, ['Q'] = 13, ['r'] = 14, ['R'] = 14, ['s'] = 15, ['S'] = 15,
  ['t'] = 16, ['T'] = 16, ['v'] = 17, ['V'] = 17, ['w'] = 18, ['W'] = 18, ['y'] = 19, ['Y'] = 19
};
const char Mapbase[] = { 'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' };

#else

const int Basemap[128] =
{
  [0 ... 127] = DEFAULT,
  ['a'] = 0, ['A'] = 0,
  ['c'] = 1, ['C'] = 1,
  ['g'] = 2, ['G'] = 2,
  ['t'] = 3, ['T'] = 3,
};
const char Mapbase[]={'A','C','G','T'};

#endif

//primers for hashtable size
const unsigned int primer[25] =
{ //hashsize nearest larger power of 2:  2^8 ~ 2^32, index 0 ~ 24
  251, 509, 1021, 2039, 4093, 8191, 16381,
  32749, 65521, 131071, 262139, 524287,
  1048573, 2097143, 4194301, 8388593, 16777213,
  33554393, 67108859, 134217689, 268435399,
  536870909, 1073741789, 2147483647, 4294967291
};

/*system status*/
double get_sys_mmry(void)
{
  struct sysinfo myinfo;
  double available_bytes;// total_bytes;
  sysinfo(&myinfo) ;
	available_bytes = myinfo.mem_unit * myinfo.totalram ;//+ myinfo.bufferram) ;
	return (available_bytes/BBILLION);
};


/*allowed input file format*/
const char *acpt_infile_fmt[ACPT_FMT_SZ] = {
	"fna",
	"fas",
	"fasta",
	"fq",
	"fastq",
	"fa",
	"co"
};
//stage I fmt
const char *fasta_fmt[FAS_FMT_SZ] = {
	"fasta",
	"fna",
	"fas",
	"fa"	
};
const char *fastq_fmt[FQ_FMT_SZ] = {
	"fq",
	"fastq"
};
// stage II format
const char *co_fmt[CO_FMT_SZ] = {
	"co"
};
// stage III format
const char *mco_fmt[MCO_FMT_SZ] = {
	"mco"
};
// compresssion format
const char *compress_fmt[CMPRESS_FMT_SZ] = {
  ".gz",
  ".bz2"
};

/** Logging **/
void log_printf(struct arg_global* g, int level, const char* fmt, ...)
{
  va_list ap;
  FILE* f = stdout;

  if(g->verbosity < level)
    return;

  if(level == 0)
    f = stderr;

  va_start(ap, fmt);

  vfprintf(f, fmt, ap);

  va_end(ap);
}

FILE * fpathopen (const char *dpath, const char *fname, const char *mode )
{
  char *fullname = malloc(PATHLEN*sizeof(char));
  sprintf(fullname,"%s/%s",dpath,fname);
	
	struct stat s;
	if(! ( ( stat(dpath, &s) == 0 ) && S_ISDIR(s.st_mode) ) )	
		mkdir(dpath, 0777);

  FILE *fp;
  if( (fp = fopen(fullname, mode) ) == NULL )
    err(errno,"fpathopen()::%s",fullname);

  return fp;
}

/*read file list into organized array********
* detect bad path(too long,empty row, not a file et al.) 
**********************/
//20190910, enable option format check for sra accession format
infile_tab_t * organize_infile_list(char* list_path, int fmt_ck) //20190910, enable option format check
{	
		infile_tab_t *infile_stat = malloc(sizeof(infile_tab_t)); 
		int alloc_usize = 1024;
		infile_stat->organized_infile_tab = malloc( sizeof(infile_entry_t) * alloc_usize );
    struct stat path_stat;
		
    FILE * list;
    list = fopen(list_path,"r");
    if(!list) err(errno,"can't open file %s",list_path);
    char *buf = malloc( LMAX * sizeof(char));
    int file_num = 0;

		if(fmt_ck){ //format check
    	while ( (fgets(buf,LMAX,list))!=NULL){
				while(isspace(*buf)) buf++; //rm spaces before path 
      	buf[strcspn(buf, "\r\n")] = 0; //trim \r\n
      	if( strlen(buf) < 1 )
        	continue;
      	if( strlen(buf) > PATHLEN )
        	err(errno,"the input list: %s\n %dth line:  %s  has %lu characters exceed the maximal allowed length %d",
            list_path,file_num, buf,strlen(buf),PATHLEN);
      //reset path_stat everytime before stat new file
      	memset(&path_stat, 0, sizeof path_stat);
      	stat(buf, &path_stat) ;
      	if(!S_ISREG(path_stat.st_mode))
          err(errno,"%dth line: %s",file_num, buf);
				else if(!isOK_fmt_infile(buf,acpt_infile_fmt,ACPT_FMT_SZ)){				
						printf ("isOK_fmt_infile(): wrong format %dth line: %s\nSupported format are:\n",file_num, buf);
						for(int i=0; acpt_infile_fmt[i]!=NULL;i++)
							printf(".%s ",acpt_infile_fmt[i]);
						printf("\n");	
							err(errno,"program exit");
				}
				else {
					infile_stat->organized_infile_tab[file_num].fsize = path_stat.st_size;
					infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char) );
					strcpy(infile_stat->organized_infile_tab[file_num].fpath, buf);
					file_num++;
					if( file_num >= alloc_usize){
      			alloc_usize+=alloc_usize;
        		infile_stat->organized_infile_tab
                = realloc(infile_stat->organized_infile_tab, sizeof(infile_entry_t) * alloc_usize );
       		}
				}
    	};// while 
		}//if 
		else{ // no fmt check
			while ( (fgets(buf,LMAX,list))!=NULL){
        while(isspace(*buf)) buf++; //rm spaces before path
        buf[strcspn(buf, "\r\n")] = 0; //trim \r\n
        if( strlen(buf) < 1 )
          continue;
        if( strlen(buf) > PATHLEN )
          err(errno,"the input list: %s\n %dth line:  %s  has %lu characters exceed the maximal allowed length %d",
            list_path,file_num, buf,strlen(buf),PATHLEN);

				infile_stat->organized_infile_tab[file_num].fsize = 0;
				infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char) );
				strcpy(infile_stat->organized_infile_tab[file_num].fpath, buf);
				file_num++;
				if( file_num >= alloc_usize){
            alloc_usize+=alloc_usize;
            infile_stat->organized_infile_tab
                = realloc(infile_stat->organized_infile_tab, sizeof(infile_entry_t) * alloc_usize );
       }
			
		}//while
	} //else


    infile_stat->infile_num = file_num ;
    fclose(list);
		free(buf);
    return infile_stat;
};
//20190910: enhanced option for format check
infile_tab_t * organize_infile_frm_arg (int num_remaining_args, char ** remaining_args,int fmt_ck)
{
//	printf("organize_infile_frm_arg:%d\t%s\n", num_remaining_args,remaining_args[0]);
	infile_tab_t *infile_stat = malloc(sizeof(infile_tab_t));
	int file_num = 0;	
	struct stat path_stat;
	DIR *dirp;
	struct dirent *dirent;
	char fullpath[PATHLEN];
	int alloc_usize = 1024;
	infile_stat->organized_infile_tab = malloc( sizeof(infile_entry_t) * alloc_usize );

if(fmt_ck) { //normal case: files
	for(int i=0;i<num_remaining_args;i++){
		stat(remaining_args[i],&path_stat);
		if( S_ISDIR(path_stat.st_mode)){
				if (( dirp = opendir(remaining_args[i]) ) == NULL)
					err(errno, "%dth argument: can't open %s",i+1, remaining_args[i] );

				while ((dirent = readdir(dirp)) != NULL){
					if(strlen(remaining_args[i]) + strlen(dirent->d_name) + 1 >PATHLEN) 
						err(errno,"path: %s/%s exceed maximal path lenth %d",remaining_args[i], dirent->d_name,PATHLEN);
					sprintf(fullpath, "%s/%s", remaining_args[i], dirent->d_name);
					stat(fullpath,&path_stat);
			
      		if(isOK_fmt_infile(fullpath,acpt_infile_fmt,ACPT_FMT_SZ)){						
						infile_stat->organized_infile_tab[file_num].fsize = path_stat.st_size;
						infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char));
						sprintf(infile_stat->organized_infile_tab[file_num].fpath, "%s/%s", remaining_args[i], dirent->d_name);
						file_num++;
						if( file_num >= alloc_usize){
							alloc_usize+=alloc_usize;

							infile_stat->organized_infile_tab 
								= realloc(infile_stat->organized_infile_tab,sizeof(infile_entry_t) * alloc_usize );
						} 						
					}
				}
				closedir(dirp);
		}
		else if(isOK_fmt_infile(remaining_args[i],acpt_infile_fmt,ACPT_FMT_SZ)){
			stat(remaining_args[i],&path_stat);		
			infile_stat->organized_infile_tab[file_num].fsize = path_stat.st_size;
			infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char));
			strcpy(infile_stat->organized_infile_tab[file_num].fpath, remaining_args[i]);
			file_num++;
			if( file_num >= alloc_usize){
      	alloc_usize+=alloc_usize;

        infile_stat->organized_infile_tab 
					= realloc( infile_stat->organized_infile_tab, sizeof(infile_entry_t) * alloc_usize );
      }
		}
		else{
			printf ("wrong format %dth argument: %s\nSupported format are:\n",i+1,remaining_args[i]);
      for(int i=0; acpt_infile_fmt[i]!=NULL;i++)
      	printf(".%s ",acpt_infile_fmt[i]);
      printf("\n");
      err(errno,"program exit");
		}
	};
}// fmt check mode end
else{ // accession case non-file mode
	for(int i=0;i<num_remaining_args;i++){
		infile_stat->organized_infile_tab[file_num].fsize = 0;
		infile_stat->organized_infile_tab[file_num].fpath = malloc(PATHLEN * sizeof(char));
		strcpy(infile_stat->organized_infile_tab[file_num].fpath, remaining_args[i]);
		file_num++;
		if( file_num >= alloc_usize){
        alloc_usize+=alloc_usize;

        infile_stat->organized_infile_tab
          = realloc( infile_stat->organized_infile_tab, sizeof(infile_entry_t) * alloc_usize );
    }
	}
}

	infile_stat->infile_num = file_num ;
	return infile_stat;
};

/*count file for each fmt type*/
infile_fmt_count_t *infile_fmt_count ( infile_tab_t * infile_tab )
{
	infile_fmt_count_t tmp_fmt_count = {0,0,0,0};
	infile_fmt_count_t* fmt_count = (infile_fmt_count_t*)malloc(sizeof(infile_fmt_count_t));
	*fmt_count = tmp_fmt_count;
	for(int i = 0; i < infile_tab->infile_num; i++ )
	{

		if( isOK_fmt_infile(infile_tab->organized_infile_tab[i].fpath,fasta_fmt,FAS_FMT_SZ) )
			fmt_count->fasta++;
		else if (isOK_fmt_infile (infile_tab->organized_infile_tab[i].fpath, fastq_fmt,FQ_FMT_SZ ) )
			fmt_count->fastq++;
		else if(isOK_fmt_infile (infile_tab->organized_infile_tab[i].fpath, co_fmt,CO_FMT_SZ) )
			fmt_count->co++;
		else if(isOK_fmt_infile (infile_tab->organized_infile_tab[i].fpath, mco_fmt,MCO_FMT_SZ))
			fmt_count->mco++;
		else if(infile_tab->organized_infile_tab[i].fsize != 0) 
			err(errno,"infile_fmt_count(): %s is not accept format(.fasta,.fastq,.co)",infile_tab->organized_infile_tab[i].fpath);
	}
	return fmt_count;
}

/*bin stat, get basename for mco binning*/
#define GZCOMPRESS_RATE 5
bin_stat_t * get_bin_basename_stat( infile_entry_t* organized_infile_tab, int *shuffle_arr,int binsz)
{
	bin_stat_t *ret = malloc(sizeof(bin_stat_t));	
	ret->seqfilebasename = malloc( binsz * BASENAME_LEN ); 		
	llong fsize;
	ret->est_kmc_bf_dr = 0;

	char *fullname,*filename;
	char suftmp[10];
	char cp_filename[BASENAME_LEN];
  int compress_fmt_count = sizeof(compress_fmt)/sizeof(compress_fmt[0]) ;
	int acpt_infile_fmt_count = sizeof(acpt_infile_fmt)/sizeof(acpt_infile_fmt[0]);
	int basename_len;

	for(int i=0; i < binsz; i++){
		fsize = organized_infile_tab[shuffle_arr[i]].fsize;
		fullname = organized_infile_tab[shuffle_arr[i]].fpath;
		(filename = strrchr(fullname,'/') ) ? ++filename : (filename = fullname);		

		if( strlen(filename) > BASENAME_LEN)
			err(errno,"input filename:%s excess %d ",filename,BASENAME_LEN);
		strcpy(cp_filename,filename);		

		// if file is compressed
		for(int j = 0; j < compress_fmt_count;j++ ){

			basename_len = strlen(filename) - strlen(compress_fmt[j]);
			if( strcmp( ( filename + basename_len ),compress_fmt[j]) == 0 ){
				strcpy( ret->seqfilebasename[i],filename);
				*(cp_filename + basename_len) = '\0'; //trunck .gz,.bz2 .. 
				fsize *= GZCOMPRESS_RATE; // estimate original file size
				break;
			};
		};

		for(int j = 0; j < acpt_infile_fmt_count; j++){
			sprintf(suftmp,".%s",acpt_infile_fmt[j]);
			basename_len = strlen(cp_filename) - strlen(suftmp);
			if( strcmp( (cp_filename + basename_len ), suftmp ) == 0 ){
				if( isOK_fmt_infile(cp_filename,fastq_fmt,FQ_FMT_SZ) )
					fsize = fsize / 2 ; //quality occupy half of the file
				else if( isOK_fmt_infile(cp_filename, co_fmt,CO_FMT_SZ ) )
					fsize = fsize / sizeof(llong) ; 
				*(cp_filename + basename_len) = '\0';
				break;
			}
		};
		//estimated kmer count sum accross files in the bin
		ret->est_kmc_bf_dr += fsize; 
		strcpy(ret->seqfilebasename[i],cp_filename); 	
	};
	
	return ret;
}

int str_suffix_match(char *str, const char *suf)
{
	int ret = 0;
	if( (strlen(str) > strlen(suf)) && (strcmp( (str + strlen(str) - strlen(suf)),suf) == 0) )
		ret = 1;
	return ret;
};

const char * get_pathname(const char *fullpath, const char *suf)
{
	char *pathcp = malloc( strlen(fullpath) + 1) ;
	strcpy(pathcp,fullpath);
	*(pathcp + strlen(pathcp) - strlen(suf)) = '\0';
	return pathcp; 
}


llong find_lgst_primer_2pow(int w)
{
	if( w < 2 || w > 62 ){
		perror("find_1st_primer_after_2pow: argument should between 8 and 62");
		exit(EXIT_FAILURE);
	}

 	llong n = ( 1llu << w ) ; //memory limit for possible subcontext space
	llong hshsz = (llong) ( (double)n * CTX_SPC_USE_L / LD_FCTR) ; 
	printf("w=%d\tspace_sz=%llu\thashsize=%llu\tkmerlimt=%llu\n",w,n,hshsz,(llong)(hshsz*LD_FCTR) ) ; 
	llong i = 3, c ; llong prime = 0;
	for(i = n - 1 ; i > (n >> 1) ; i--) 
	{
		for ( c = 2  ; c <= (int)pow(i+1,0.5) ; c++ )
		{
			if( i%c == 0 )
				break;
		}
		
		if( c*c > i ){
			prime = i;
			break;
		}
	} 
	printf("nearest prime=%llu\n",prime);
	return prime;
}


// Helper function to check if a number is prime
static uint32_t is_prime(uint32_t m) {
    if (m <= 1) return 0;
    if (m == 2) return 1;
    if (m % 2 == 0) return 0;
    
    uint32_t sqrt_m = (uint32_t)sqrt(m);
    for (uint32_t i = 3; i <= sqrt_m; i += 2) {
        if (m % i == 0) return 0;
    }
    return 1;
}

// Main function to find the next prime after n
uint32_t nextPrime(uint32_t n) {
    if (n < 2) return 2;
    if (n == 2) return 3;
    
    // Start checking from the next odd number
    uint32_t candidate = (n % 2 == 0) ? n + 1 : n + 2;
    
    while (1) {
        if (is_prime(candidate)) {
            return candidate;
        }
        candidate += 2; // Skip even numbers
    }
}

void replaceChar(char *str, char oldChar, char newChar) {
    while (*str != '\0') {
        if (*str == oldChar) {
            *str = newChar;
        }
        str++;
    }
}


#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>

#ifdef _WIN32
#include <direct.h>
#define MKDIR(path) _mkdir(path)
#else
#include <unistd.h>
#define MKDIR(path) mkdir(path, 0755)
#endif

int mkdir_p(const char *path) {
	char temp[1024];
    size_t len = strlen(path);
    struct stat st;

    // Copy the path to a temporary buffer
    snprintf(temp, sizeof(temp), "%s", path);

    // Remove trailing slash if present
    if (temp[len - 1] == '/')
        temp[len - 1] = '\0';

    // Create directories recursively
    for (char *p = temp + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0'; // Temporarily terminate the string
            if (stat(temp, &st) != 0) { // If the directory doesn't exist
                if (MKDIR(temp) != 0 && errno != EEXIST) {
                    perror("mkdir");
                    return -1;
                }
            }
            *p = '/'; // Restore the slash
        }
    }

    // Create the final directory
    if (stat(temp, &st) != 0) {
        if (MKDIR(temp) != 0 && errno != EEXIST) {
            perror("mkdir");
            return -1;
        }
    }
    return 0; // Success
}

//test if parent_path contain file dstat_f, return full path if has otherwise return NULL
//orginally develop fo test dstat_f , but can generally test other file
char* test_get_fullpath(const char *parent_path, const char *dstat_f)
{
  struct stat path_stat;
  if( stat(parent_path, &path_stat) < 0 )
    err(errno,"test_get_fullpath()::%s",parent_path);
  if( S_ISDIR(path_stat.st_mode) ){
    char* fullpath = malloc(PATHLEN+1);
    sprintf((char *)fullpath,"%s/%s", parent_path, dstat_f);
    FILE *fp;
    if ( (fp = fopen(fullpath,"rb")) != NULL ){
      fclose(fp);
      return fullpath;
    }
    else{
	  printf("%s(): %s do not exists\n",__func__,fullpath);
      free((char*)fullpath);
      return NULL;
    }
  }
  else{
		printf("%s()::%s is not a director\n",__func__,parent_path);
    return NULL;
	}
};

char* test_create_fullpath(const char *parent_path, const char *dstat_f){
	mkdir_p(parent_path);
	char* fullpath = malloc(PATHLEN+1);
	sprintf((char *)fullpath,"%s/%s", parent_path, dstat_f);

	 FILE *fp;
   if ( (fp = fopen(fullpath,"wb")) != NULL ){
      fclose(fp);
      return fullpath;
    }
    else{
      free((char*)fullpath);
      return NULL;
    }		
}

int file_exists_in_folder(const char *folder, const char *filename) {
    char filepath[1024]; // Buffer to store the full file path
    snprintf(filepath, sizeof(filepath), "%s/%s", folder, filename); // Construct path    
    // Check if the file exists
    return access(filepath, F_OK) == 0;
}

//by chatgpt
// fread wrapper to read the entire file into memory
#define MAX_FREAD_SIZE (16LU*1024*1024*1024)
void *read_from_file(const char *file_path, size_t *file_size) {
    // Open the file in binary read mode
    FILE *file = fopen(file_path, "rb");
    if (!file) err(EXIT_FAILURE, "%s(): Error opening file '%s'", __func__, file_path);
     // Seek to the end of the file to determine its size
    if (fseek(file, 0, SEEK_END) != 0)  err(EXIT_FAILURE, "%s(): Failed to seek to the end of file '%s'", __func__, file_path); 
    // Get the file size
    long size = ftell(file);
    if (size == -1)  err(EXIT_FAILURE, "%s(): Failed to determine file size for '%s'", __func__, file_path); 
   // Seek back to the beginning of the file
    rewind(file);

    void *buffer;	
    if(size < MAX_FREAD_SIZE) {
    // Allocate memory to hold the file contents
    	if ( (buffer = malloc(size)) == NULL) {
        	fclose(file);
       		 err(EXIT_FAILURE, "%s(): Memory allocation failed for file '%s'", __func__, file_path);
    	}
    	// Read the entire file into memory
    	size_t read_count = fread(buffer, 1, size, file);
    	if (read_count != (size_t)size) {
        	free(buffer);
        	fclose(file);
        	err(EXIT_FAILURE, "%s(): Failed to read the file '%s' into memory", __func__, file_path);
    	}
    	fclose(file);
    }
    else{ //mmap
    	int fd = open(file_path, O_RDWR, 0666);
	if (fd == -1) err(EXIT_FAILURE,"%s(): open %s failed",__func__, file_path); 
	if ((buffer=mmap(NULL, size, PROT_READ, MAP_SHARED, fd, 0)) == MAP_FAILED) {
    		close(fd);
    		err(EXIT_FAILURE, "%s(): Failed to mmap the file '%s'", __func__, file_path);
	} 
        close(fd);	
    }
    // Set the file size if a valid pointer is provided
    if (file_size)  *file_size = (size_t)size;
    return buffer;
}

void free_read_from_file (void *buffer, size_t file_size){
	if(file_size < MAX_FREAD_SIZE) free(buffer);
	else if(munmap(buffer,file_size) == -1) 
		err(EXIT_FAILURE, "%s(): Failed to munmap buffer with size of '%lu'", __func__, file_size);		
}

// Function to read a file into memory (supports both text and binary modes)
void *read_file_mode(const char *file_path, size_t *file_size, const char *mode) {
    // Validate mode input
    if (!mode || (strcmp(mode, "r") != 0 && strcmp(mode, "rb") != 0)) 
        err(EXIT_FAILURE, "%s(): Invalid mode '%s'. Use \"r\" for text or \"rb\" for binary.", __func__, mode ? mode : "(null)");
    // Open the file with the specified mode
    FILE *file = fopen(file_path, mode);
    if (!file)  err(EXIT_FAILURE, "%s(): Error opening file '%s': %s", __func__, file_path, strerror(errno));
     // Seek to the end of the file to determine its size
    if (fseek(file, 0, SEEK_END) != 0) {
		fclose(file);
        err(EXIT_FAILURE, "%s(): Failed to seek to the end of file '%s'", __func__, file_path);
    }
    // Get the file size
    long size = ftell(file);
    if (size == -1) {
		fclose(file);
        err(EXIT_FAILURE, "%s(): Failed to determine file size for '%s'", __func__, file_path);
    }
    // Seek back to the beginning of the file
    rewind(file);
    // Check if text mode ("r") needs an extra null-terminator
    int is_text = (strcmp(mode, "r") == 0);
    void *buffer = malloc(size + (is_text ? 1 : 0));
    if (!buffer) {
		fclose(file);
        err(EXIT_FAILURE, "%s(): Memory allocation failed for file '%s'\n", __func__, file_path);
    }
    // Read the entire file into memory
    size_t read_count = fread(buffer, 1, size, file);
    if (read_count != (size_t)size) {
        free(buffer);
        fclose(file);
        err(EXIT_FAILURE, "%s(): Failed to read the file '%s' into memory\n", __func__, file_path); 
    }

    fclose(file);
    // Null-terminate if reading as a text file
    if (is_text) ((char *)buffer)[size] = '\0';
    // Set the file size if a valid pointer is provided
    if (file_size) *file_size = (size_t)size;

    return buffer;
}
// Function to read a file and split it into an array of lines
char **read_lines_from_file(const char *file_path, int *line_count) {
    size_t file_size;
    char *content = (char *)read_file_mode(file_path, &file_size, "r");  // Read file as text
    if (!content) return NULL;
    // Count the number of lines
    *line_count = 0;
    for (size_t i = 0; i < file_size; i++) {
        if (content[i] == '\n') (*line_count)++;
    }
    (*line_count)++;  // Last line (if file doesn't end with '\n')
    // Allocate an array to hold pointers to each line
    char **lines = malloc((*line_count) * sizeof(char *));
    if (!lines) {
        free(content);
        err(EXIT_FAILURE, "%s():Memory allocation failed for lines array",__func__);
    }
    // Tokenize the file content into lines
    size_t index = 0;
    char *line = strtok(content, "\n");
    while (line) {
        lines[index++] = line;
        line = strtok(NULL, "\n");
    }
    return lines;  // Note: The returned array points to `content` memory
}

// fwrite wrapper to write memory to a file
void write_to_file(const char *file_path, const void *data, size_t data_size) {
    // Open the file in binary write mode
    FILE *file = fopen(file_path, "wb");
    if (!file) {
        err(EXIT_FAILURE, "%s(): Error opening file '%s' for writing", __func__, file_path);
    }
    // Write the data to the file
    size_t write_count = fwrite(data, 1, data_size, file);
    if (write_count != data_size) {
        fclose(file);
        err(EXIT_FAILURE, "%s(): Failed to write data to file '%s'", __func__, file_path);
    }

    // Close the file
    if (fclose(file) != 0) {
        err(EXIT_FAILURE, "%s(): Failed to close file '%s'", __func__, file_path);
    }
}

void concat_and_write_to_file(const char *file_path, const void *block1, size_t size1, const void *block2, size_t size2) {
    // Calculate the total size of the new block
    size_t total_size = size1 + size2;

    // Allocate memory for the new block
    void *joint_block = malloc(total_size);
    if (!joint_block) {
        err(EXIT_FAILURE, "concat_and_write_to_file(): Memory allocation failed");
    }

    // Copy the first block into the new block
    memcpy(joint_block, block1, size1);

    // Append the second block after the first
    memcpy((char *)joint_block + size1, block2, size2);

    // Write the concatenated block to the file
    write_to_file(file_path, joint_block, total_size);

    // Free the allocated memory
    free(joint_block);
}


//vector

// Initialize the vector
void vector_init(Vector *vec, size_t element_size) {
    vec->data = NULL;
    vec->element_size = element_size;
    vec->size = 0;
    vec->capacity = 0;
}

// Free the vector's memory
void vector_free(Vector *vec) {
    free(vec->data);
    vec->data = NULL;
    vec->size = 0;
    vec->capacity = 0;
}

// Add an element to the vector
void vector_push(Vector *vec, const void *element) {
    if (vec->size == vec->capacity) {
        size_t new_capacity = vec->capacity == 0 ? 4 : vec->capacity * 2;
        vec->data = realloc(vec->data, new_capacity * vec->element_size);
        if (vec->data == NULL) {
            fprintf(stderr, "Failed to allocate memory\n");
            exit(EXIT_FAILURE);
        }
        vec->capacity = new_capacity;
    }
    memcpy((char *)vec->data + vec->size * vec->element_size, element, vec->element_size);
    vec->size++;
}

// Get an element from the vector (by index)
void *vector_get(Vector *vec, size_t index) {
    if (index >= vec->size) {
        fprintf(stderr, "Index out of bounds\n");
        exit(EXIT_FAILURE);
    }
    return (char *)vec->data + index * vec->element_size;
} 
void vector_reserve(Vector *v, size_t new_capacity) {
    if (new_capacity > v->capacity) {
        void *new_data = realloc(v->data, new_capacity * v->element_size);
        if (!new_data) exit(EXIT_FAILURE);
        v->data = new_data;
        v->capacity = new_capacity;
    }
}
//free all
#include <stdarg.h>
void free_all(void *first, ...) {
    va_list args;
    va_start(args, first);

    void *ptr = first;
    while (ptr != NULL) {
        free(ptr);
        ptr = va_arg(args, void *);
    }

    va_end(args);
}

// Wrapper for snprintf to return the concatenated string directly
char* format_string(const char* format, ...) {
    va_list args;

    // Step 1: Determine the length of the formatted string
    va_start(args, format);
    int len = vsnprintf(NULL, 0, format, args); // Get required buffer size
    va_end(args);

    if (len < 0) {
        return NULL; // vsnprintf failed
    }

    // Step 2: Allocate memory for the formatted string
    char* result = (char*)malloc((len + 1) * sizeof(char)); // +1 for null terminator
    if (!result) {
        return NULL; // Memory allocation failed
    }

    // Step 3: Format the string into the allocated buffer
    va_start(args, format);
    vsnprintf(result, len + 1, format, args); // Write the formatted string
    va_end(args);

    return result;
}

unify_sketch_t* generic_sketch_parse(const char *qrydir) {
    unify_sketch_t *result = calloc(1,sizeof(unify_sketch_t)); // Dynamically allocate memory
    if (result == NULL) {
        err(EXIT_FAILURE, "%s(): Memory allocation failed", __func__);
    }
    size_t file_size, combco_fsize ;

    if (file_exists_in_folder(qrydir, co_dstat)) {
        result->stat_type = 1;
        result->mem_stat = read_from_file(test_get_fullpath(qrydir, co_dstat), &file_size);
        memcpy(&result->stats.co_stat_val, result->mem_stat, sizeof(co_dstat_t));
        if (result->stats.co_stat_val.comp_num > 1) {
            err(EXIT_FAILURE, "%s(): comp_num must be 1. Found: %d", __func__, result->stats.co_stat_val.comp_num);
        }

        result->hash_id = result->stats.co_stat_val.shuf_id;
        result->infile_num = result->stats.co_stat_val.infile_num;
        result->kmerlen = result->stats.co_stat_val.kmerlen;
        result->gname = (char (*)[PATHLEN])(result->mem_stat + sizeof(co_dstat_t) + sizeof(ctx_obj_ct_t) * result->infile_num);

        unsigned int *tmp_combco = read_from_file(format_string("%s/%s.0", qrydir, skch_prefix), &combco_fsize);
  	size_t *tmp_index_combco = read_from_file(format_string("%s/%s.0", qrydir, idx_prefix), &file_size);
        result->comb_sketch = malloc(sizeof(uint64_t) * tmp_index_combco[result->infile_num]);
        for (uint64_t i = 0; i < tmp_index_combco[result->infile_num]; i++) result->comb_sketch[i] = tmp_combco[i];
		free_read_from_file(tmp_combco,combco_fsize);
        
		if(result->stats.co_stat_val.koc) {
			unsigned short *tmp_ab = read_from_file(format_string("%s/%s.0.a", qrydir, skch_prefix), &file_size);
			result->abundance =  malloc(sizeof(uint32_t) * tmp_index_combco[result->infile_num]);
			for (uint64_t i = 0; i < tmp_index_combco[result->infile_num]; i++) result->abundance[i] = tmp_ab[i]; 
			free(tmp_ab);
		}

        if (sizeof(size_t) == sizeof(uint64_t)) {
            result->sketch_index = (uint64_t *)tmp_index_combco;
        } else {
            result->sketch_index = malloc(sizeof(uint64_t) * (result->infile_num + 1));
            for (int i = 0; i < result->infile_num + 1; i++) {
                result->sketch_index[i] = tmp_index_combco[i];
            }
            free(tmp_index_combco);
        }
    } else if (file_exists_in_folder(qrydir, sketch_stat)) {
        result->stat_type = 2;
        result->mem_stat = read_from_file(test_get_fullpath(qrydir, sketch_stat), &file_size);
        memcpy(&result->stats.lco_stat_val, result->mem_stat, sizeof(dim_sketch_stat_t));

        result->hash_id = result->stats.lco_stat_val.hash_id;
        result->infile_num = result->stats.lco_stat_val.infile_num;
        result->kmerlen = result->stats.lco_stat_val.klen;
        result->gname = (char (*)[PATHLEN])(result->mem_stat + sizeof(dim_sketch_stat_t));
        result->comb_sketch = read_from_file(test_get_fullpath(qrydir, combined_sketch_suffix), &file_size);
        result->sketch_index = read_from_file(test_get_fullpath(qrydir, idx_sketch_suffix), &file_size);
		if(result->stats.lco_stat_val.koc) result->abundance = read_from_file(test_get_fullpath(qrydir, combined_ab_suffix), &file_size);

    } else {
        err(EXIT_FAILURE, "%s(): No valid sketch data found in folder: %s", __func__, qrydir);
    }

    return result; // Return the dynamically allocated memory
}

void free_unify_sketch (unify_sketch_t *result) {
    if (result == NULL) return;  // Avoid dereferencing a NULL pointer   
    free_read_from_file(result->comb_sketch,sizeof(result->comb_sketch[0])*result->sketch_index[result->infile_num] );
    free_all(result->sketch_index,result->mem_stat,NULL);
    result->mem_stat = result->sketch_index = result->comb_sketch = NULL;
	if(result->abundance) free(result->abundance);
	free(result);
}


void replace_special_chars_with_underscore(char *str) {
    for (int i = 0; str[i] != '\0'; i++) {
        if (!isalnum(str[i])) {  // If the character is not alphanumeric
            str[i] = '_';         // Replace it with an underscore
        }
    }
}


uint64_t GetAvailableMemory() {
    struct sysinfo info;
    if (sysinfo(&info) == -1) {
        perror("sysinfo failed");
        return 0;
    }
    
    return (info.freeram * info.mem_unit);
}
