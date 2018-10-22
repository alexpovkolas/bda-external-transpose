#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>


#define __PROFILE__

#ifdef __PROFILE__


#endif

using namespace std;

void transpose(char *dst, const char *src, size_t n, size_t m) {
    size_t block = 32;
    for (size_t i = 0; i < n; i += block) {
        for(size_t j = 0; j < m; ++j) {
            for(size_t b = 0; b < block && i + b < n; ++b) {
                dst[j*n + i + b] = src[(i + b)*m + j];
            }
        }
    }
}

void external_transpose(ifstream &in, ofstream &out, int n, int m, int memory_size) {
    int read_block_size = memory_size / 2 - 1000;

    int lines_block_count = 0;    // lines count we read for one block
    int cols_block_count = sqrt(read_block_size); // columns count we read for one block. One block is matrix lines_count * cols_count
    if (cols_block_count < n && cols_block_count < m) {
        lines_block_count = cols_block_count;
    } else if (cols_block_count > n && cols_block_count > m){
        lines_block_count = n;
        cols_block_count = m;
    } else if (cols_block_count > n) {
        lines_block_count = n;
        cols_block_count = read_block_size / lines_block_count;
    } else if (cols_block_count > m) {
        cols_block_count = m;
        lines_block_count = read_block_size / cols_block_count;
    }


    char *source = new char[read_block_size];
    char *destination = new char[read_block_size];
    int params_offset = 8;

    for (int i = 0; i < ceil(n / (double)lines_block_count); ++i) {
        for (int j = 0; j < ceil(m / (double)cols_block_count); ++i) {


            // read block
            int lines_to_read = n - lines_block_count * i > lines_block_count ? lines_block_count : n - lines_block_count * i;
            int cols_to_read = m - cols_block_count * j > cols_block_count ? cols_block_count : n - cols_block_count * i;

            for (int k = 0; k < lines_to_read; ++k) {
                long source_offset = params_offset +         // n and m 8 bytes
                        (k + i * lines_block_count) * m +    // lines offset
                                     j * cols_block_count;   // offset inside current line
                in.seekg(source_offset);
                in.read(source + k * cols_block_count, cols_to_read);
            }

            transpose(destination, source, lines_to_read, cols_to_read);


            // write block
            for (int k = 0; k < cols_to_read; ++k) {
                // Mirror source offset
                long dest_offset =  (k + j * cols_block_count) * n +
                                     i * lines_block_count;
                out.seekp(dest_offset);
                out.write(source + k * lines_block_count, lines_to_read);
            }
        }
    }
}

int main() {

    ifstream in("input.bin", ios::binary);
    ofstream out("output.bin", ios::binary | ios::out);

    int n = 0;
    int m = 0;

    in >> n >> m;

    int memory_limit = 1000000;
    external_transpose(in, out, n, m, memory_limit);

    return 0;
}
