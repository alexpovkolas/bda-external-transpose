#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>


//#define __PROFILE__

#ifdef __PROFILE__

#include <algorithm>
#include <chrono>

#endif


using namespace std;

void transpose(char *dst, const char *src, size_t n, size_t m) {
    size_t block = 16;
    for (size_t i = 0; i < n; i += block) {
        for(size_t j = 0; j < m; ++j) {
            for(size_t b = 0; b < block && i + b < n; ++b) {
                dst[j*n + i + b] = src[(i + b)*m + j];
            }
        }
    }
}

void external_transpose(ifstream &in, ofstream &out, int n, int m, int memory_size) {
    int read_block_size = memory_size / 2;

    int lines_block_count = 0;    // lines count we read for one block
    int cols_block_count = sqrt(read_block_size); // columns count we read for one block. One block is matrix lines_count * cols_count
    if (cols_block_count < n && cols_block_count < m) {
        lines_block_count = cols_block_count;
    } else if (cols_block_count >= n && cols_block_count >= m){
        lines_block_count = n;
        cols_block_count = m;
    } else if (cols_block_count >= n) {
        lines_block_count = n;
        cols_block_count = read_block_size / lines_block_count;
    } else if (cols_block_count >= m) {
        cols_block_count = m;
        lines_block_count = read_block_size / cols_block_count;
    }

#ifdef __PROFILE__

    int file_reads_count = 0;
    int file_writes_count = 0;

    cout << "lines_block_count = " << lines_block_count << " cols_block_count = " << cols_block_count << endl;
#endif



    char *source = new char[lines_block_count * cols_block_count];
    char *destination = new char[lines_block_count * cols_block_count];
    int params_offset = 8;

    int vert_blocks = ceil(n / (double)lines_block_count);
    int hor_blocks = ceil(m / (double)cols_block_count);

#ifdef __PROFILE__
    cout << "vert_blocks = " << vert_blocks << " hor_blocks = " << hor_blocks << endl;
#endif

    for (int i = 0; i < vert_blocks; ++i) {
        for (int j = 0; j < hor_blocks; ++j) {

#ifdef __PROFILE__
            cout << "i = " << i << " j = " << j << endl;
#endif

            // read block
            int lines_to_read = n - lines_block_count * i > lines_block_count ? lines_block_count : n - lines_block_count * i;
            int cols_to_read = m - cols_block_count * j > cols_block_count ? cols_block_count : m - cols_block_count * j;

#ifdef __PROFILE__
            cout << "lines_to_read = " << lines_to_read << " cols_to_read = " << cols_to_read << endl;
#endif

            // we can get few lines for one read operation
            if (hor_blocks == 1) {
                long source_offset = params_offset +         // n and m 8 bytes
                                     (i * lines_block_count) * m +    // lines offset
                                     j * cols_block_count;   // offset inside current line
                in.seekg(source_offset);
                in.read(source, cols_to_read * lines_to_read);

#ifdef __PROFILE__
                file_reads_count++;
                cout << "source_offset = " << source_offset << " read bytes = " << cols_to_read << endl;
#endif
            } else {
                for (int k = 0; k < lines_to_read; ++k) {
                    long source_offset = params_offset +         // n and m 8 bytes
                                         (k + i * lines_block_count) * m +    // lines offset
                                         j * cols_block_count;   // offset inside current line
                    in.seekg(source_offset);
                    in.read(source + k * cols_to_read, cols_to_read);

#ifdef __PROFILE__
                    file_reads_count++;
                    cout << "source_offset = " << source_offset << " read bytes = " << cols_to_read << endl;
#endif

                }
            }


            transpose(destination, source, lines_to_read, cols_to_read);

            // we can set few lines for one write operation
            if (vert_blocks == 1) {
                long dest_offset = params_offset +
                                   (j * cols_block_count) * n +
                                   i * lines_block_count;
                out.seekp(dest_offset);
                out.write(destination, lines_to_read * cols_to_read);

#ifdef __PROFILE__
                file_writes_count++;
                cout << "dest_offset = " << dest_offset << " write bytes = " << lines_to_read << endl;
#endif
            } else {
                // write block
                for (int k = 0; k < cols_to_read; ++k) {
                    // Mirror source offset
                    long dest_offset = params_offset +
                                       (k + j * cols_block_count) * n +
                                       i * lines_block_count;
                    out.seekp(dest_offset);
                    out.write(destination + k * lines_to_read, lines_to_read);

#ifdef __PROFILE__
                    file_writes_count++;
                    cout << "dest_offset = " << dest_offset << " write bytes = " << lines_to_read << endl;
#endif
                }
            }
        }
    }

#ifdef __PROFILE__
    cout << "file_writes_count = " << file_writes_count << " file_reads_count = " << file_reads_count << endl;
#endif

}

#ifdef __PROFILE__

void gen_test(int n, int m) {

    ofstream file("input.bin", ios::binary);
    file.write((char *)&n, 4);
    file.write((char *)&m, 4);

    for (int i = 1; i <= n * m; i++)
    {
        file.write((char *)&i, 1);
    }

    file.close();
}

bool compare_files(const string& filename1, const string& filename2)
{
    ifstream file1(filename1, ifstream::ate | ifstream::binary);
    std::ifstream file2(filename2, std::ifstream::ate | std::ifstream::binary);
    const ifstream::pos_type fileSize = file1.tellg();

    if (fileSize != file2.tellg()) {
        return false; //different file size
    }

    file1.seekg(0);
    file2.seekg(0);

    istreambuf_iterator<char> begin1(file1);
    istreambuf_iterator<char> begin2(file2);

    return equal(begin1, istreambuf_iterator<char>(),begin2); //Second argument is end-of-range iterator
}

#endif

int main() {

#ifdef __PROFILE__
//    gen_test(100000, 100); // file_writes_count = 4100 file_reads_count = 100000
//    gen_test(100, 100000); //file_writes_count = 100000 file_reads_count = 4100
//    gen_test(10000, 1000); // file_writes_count = 21000 file_reads_count = 30000
//    gen_test(1000000, 10); // file_writes_count = 410 file_reads_count = 1000000
//    gen_test(2, 10); // file_writes_count = 410 file_reads_count = 1000000
#endif


#ifdef __PROFILE__
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
#endif

    ifstream in("input.bin", ios::binary);
    ofstream out("output.bin", ios::binary | ios::out);

    int n = 0;
    int m = 0;

    in.read((char *)&n, 4);
    in.read((char *)&m, 4);
    out.write((char *)&m, 4);
    out.write((char *)&n, 4);

    int memory_limit = 500000 - 1000;

    external_transpose(in, out, n, m, memory_limit);
    out.close();

#ifdef __PROFILE__
    chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    cout << " Time difference = " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << endl;
#endif


#ifdef __PROFILE__

    cout << "Result correct - " << compare_files("output.bin", "answer.txt");

#endif


    return 0;
}
