#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
void write_matrix(const char *filename, uint16_t *data, size_t elements) {
    int fd = open(filename, O_RDWR | O_CREAT | O_TRUNC, 0644);
    size_t filesize = elements * sizeof(uint16_t);
    
    // Extend file to required size
    ftruncate(fd, filesize);
    
    // Map file to memory
    uint16_t *map = mmap(NULL, filesize, PROT_WRITE, MAP_SHARED, fd, 0);
    
    // Copy data directly to mapped memory
    memcpy(map, data, filesize);
    
    // Cleanup
    msync(map, filesize, MS_SYNC);
    munmap(map, filesize);
    close(fd);
}

/*
1.Memory-Mapped File (Fastest)
Usage:
size_t matrix_elements = (infile_num * (infile_num + 1)) / 2;
write_matrix("ctx.bin", ctx, matrix_elements);
write_matrix("obj.bin", obj, matrix_elements);
*/

#include <stdio.h>
#include <stdlib.h>

void write_matrix_blocked(FILE *fp, uint16_t *data, size_t elements) {
    const size_t block_size = 1UL << 26; // 64MB blocks
    size_t written = 0;
    
    while(written < elements) {
        size_t to_write = (elements - written > block_size) ? 
                         block_size : elements - written;
                         
        size_t ret = fwrite(data + written, sizeof(uint16_t), to_write, fp);
        if(ret != to_write) {
            perror("Write failed");
            exit(EXIT_FAILURE);
        }
        written += ret;
    }
}

/*
2. Blocked Buffered Writes (Portable)
 Usage:
FILE *f_ctx = fopen("ctx.bin", "wb");
FILE *f_obj = fopen("obj.bin", "wb");
setvbuf(f_ctx, NULL, _IOFBF, 1 << 26); // 64MB buffer
setvbuf(f_obj, NULL, _IOFBF, 1 << 26);

write_matrix_blocked(f_ctx, ctx, matrix_elements);
write_matrix_blocked(f_obj, obj, matrix_elements);

fclose(f_ctx);
fclose(f_obj);
*/


// need manually install libzstd
// make with -lzstd
/*
#include <zstd.h>

void write_compressed(const char *filename, uint16_t *data, size_t elements) {
    size_t comp_size = ZSTD_compressBound(elements * sizeof(uint16_t));
    void *compressed = malloc(comp_size);
    
    size_t result = ZSTD_compress(compressed, comp_size,
                                 data, elements * sizeof(uint16_t),
                                 ZSTD_CLEVEL_DEFAULT);
    
    FILE *fp = fopen(filename, "wb");
    fwrite(compressed, 1, result, fp);
    fclose(fp);
    free(compressed);
}
*//*
3. Compressed Storage (Space Optimized) // need manually install libzstd 
 Usage (with error checking):
write_compressed("ctx.zst", ctx, matrix_elements);
write_compressed("obj.zst", obj, matrix_elements);
*/





