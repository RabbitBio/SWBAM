#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <zlib.h>
#include <htslib/khash.h>
#include <htslib/thread_pool.h>
#include <stdint.h>
#include <chrono>
#include "config.h"
#include "BamBlock.h"
#include "Buffer.h"
#include "BamStatus.h"
#include "Duplicate.h"
#include "Overrepresent.h"
#include "CLI/CLI.hpp"
#include <sched.h>
#include <unistd.h>
#include "BamRead.h"
#include "BamCompress.h"
#include "BamCompleteBlock.h"
#include "BamTools.h"
#include "BamWriteCompress.h"
#include "BamWriter.h"
#include "BamReader.h"

#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8


#define DEBUG 0


typedef std::chrono::high_resolution_clock Clock;
#ifdef CYCLING
#define TDEF(x_) static unsigned long long int x_##_t0, x_##_t1;
#define TSTART(x_) x_##_t0 = __rdtsc();
#define TEND(x_) x_##_t1 = __rdtsc();
#define TPRINT(x_, str) printf("%-20s \t%.6f\t M cycles\n", str, (double)(x_##_t1 - x_##_t0)/1e6);
#elif defined TIMING
#define TDEF(x_) chrono::high_resolution_clock::time_point x_##_t0, x_##_t1;
#define TSTART(x_) x_##_t0 = Clock::now();
#define TEND(x_) x_##_t1 = Clock::now();
#define TPRINT(x_, str) printf("%-20s \t%.6f\t sec\n", str, chrono::duration_cast<chrono::microseconds>(x_##_t1 - x_##_t0).count()/1e6);
#else
#define TDEF(x_)
#define TSTART(x_)
#define TEND(x_)
#define TPRINT(x_, str)
#endif
//const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
//int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
//uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
//uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};


#include <iostream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>


const size_t PACK_SIZE = 1024;

const size_t QUEUE_SIZE = 1024 * 1024 * 1024;

struct bam_pack {
    std::vector<bam1_t *> bams;

    bam_pack() {
        bams.reserve(PACK_SIZE);
    }

    ~bam_pack() {
        for (auto &bam: bams) {
            bam_destroy1(bam);
        }
    }
};

std::queue<bam_pack *> bam_queue;
std::mutex queue_mutex;
std::condition_variable not_full;
std::condition_variable not_empty;

void read_thread_task(BamReader *reader) {
    long num = 0;
    bam_pack *pack = new bam_pack();
    for (size_t i = 0; i < PACK_SIZE; ++i) {
        bam1_t *b = bam_init1();
        if (b == NULL) {
            std::cerr << "[E::" << __func__ << "] Out of memory allocating BAM struct.\n";
            return;
        }
        pack->bams.push_back(b);
    }
    while (true) {
        int cnt = 0;
        for (size_t i = 0; i < PACK_SIZE; ++i) {
            if (!reader->getBam1_t(pack->bams[i])) {
                break;
            }
            num++;
            cnt++;
        }

        if (cnt == 0) {
            //delete pack;
            break;
        }

        std::unique_lock <std::mutex> lock(queue_mutex);
        not_full.wait(lock, [] { return bam_queue.size() < QUEUE_SIZE; });

        bam_queue.push(pack);
        not_empty.notify_one();
    }
    printf("num is %lld\n", num);

    std::unique_lock <std::mutex> lock(queue_mutex);
    bam_queue.push(nullptr);
    not_empty.notify_one();
}

void read_thread_task2(BamReader *reader) {
    long num = 0;
    while (true) {
        bam_pack *pack = new bam_pack();
        for (size_t i = 0; i < PACK_SIZE; ++i) {
            bam1_t *b = bam_init1();
            if (b == NULL) {
                std::cerr << "[E::" << __func__ << "] Out of memory allocating BAM struct.\n";
                delete pack;
                return;
            }

            if (!reader->getBam1_t(b)) {
                bam_destroy1(b);
                delete pack;
                break;
            }
            pack->bams.push_back(b);
        }

        if (pack->bams.empty()) {
            delete pack;
            break;
        }

        std::unique_lock <std::mutex> lock(queue_mutex);
        not_full.wait(lock, [] { return bam_queue.size() < PACK_SIZE; });

        bam_queue.push(pack);
        not_empty.notify_one();
    }

    std::unique_lock <std::mutex> lock(queue_mutex);
    bam_queue.push(nullptr);
    not_empty.notify_one();
    printf("num is %d\n", num);
}

void write_thread_task(BamWriter *writer) {
    printf("queue size %lld\n", bam_queue.size());
    double t0 = GetTime();
    while (true) {
        std::unique_lock <std::mutex> lock(queue_mutex);
        not_empty.wait(lock, [] { return !bam_queue.empty(); });

        bam_pack *pack = bam_queue.front();
        bam_queue.pop();
        //printf("queue size %lld\n", bam_queue.size());

        not_full.notify_one();

        if (pack == nullptr) {
            break;
        }

        for (auto &b: pack->bams) {
            writer->write(b);
        }
        //delete pack;
    }
    printf("main write part cost %lf\n", GetTime() - t0);
    writer->over();
}

void write_thread_task2(BamWriter *writer) {
    double t0 = GetTime();
    while (true) {
        std::unique_lock <std::mutex> lock(queue_mutex);
        not_empty.wait(lock, [] { return !bam_queue.empty(); });

        bam_pack *pack = bam_queue.front();
        bam_queue.pop();

        not_full.notify_one();

        if (pack == nullptr) {
            break;
        }

        for (auto &b: pack->bams) {
            writer->write(b);
        }
        delete pack;
    }
    printf("main write part cost %lf\n", GetTime() - t0);
    writer->over();
}

int InputBlockNum = 0;

long long NUM_N[100] = {0};
long long NUM_M[100] = {0};
long long NUM_TID[100][1000] = {0};

void read_pack_main(BGZF *fp, BamRead *read) {
    bam_block *b;
    b = read->getEmpty();
    int count = 0;
    while (read_block(fp, b) == 0) {
        read->inputBlock(b);
        b = read->getEmpty();
    }
    read->ReadComplete();
}

void write_pack(Buffer *buffer) {
    while (!buffer->is_complete()) {
        std::this_thread::sleep_for(chrono::milliseconds(10));
        buffer->output();
    }
}

void compress_pack_main(BamRead *read, BamCompress *compress) {
    pair < bam_block * , int > comp;
    bam_block *un_comp = compress->getEmpty();
    while (1) {
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
        comp = read->getReadBlock();
        //printf("%d is get One compressed data\n",id);
        if (comp.second < 0) {
            //printf("%d is Over\n",id);
            break;
        }
        block_decode_func(comp.first, un_comp);
        read->backBlock(comp.first);

        std::pair<int, int> tmp_pair = find_divide_pos_and_get_read_number(un_comp);
        un_comp->split_pos = tmp_pair.first, un_comp->bam_number = tmp_pair.second;
        compress->inputUnCompressData(un_comp, comp.second);
//        while (!compress->tryinputUnCompressData(un_comp,comp.second)){
//            std::this_thread::sleep_for(std::chrono::milliseconds(1));
//        }
        un_comp = compress->getEmpty();
    }
    compress->CompressThreadComplete();
}

void benchmark_pack(BamCompress *compress, BamCompleteBlock *completeBlock) {
    bam_block *un_comp = nullptr;
    bam_complete_block *assign_block = completeBlock->getEmpty();
    int need_block_len = 0, input_length = 0;
    int last_use_block_length = 0;
    bool isclean = true;
    int ret = -1;
    int bam_number = 0;
    while (1) {
        if (isclean && un_comp != nullptr) {
            compress->backEmpty(un_comp);
        }
        un_comp = compress->getUnCompressData();
        if (un_comp == nullptr) {
            break;
        }
        ret = un_comp->split_pos;
        need_block_len = ret;
        if (assign_block->length + need_block_len > assign_block->data_size) {
            completeBlock->backEmpty(assign_block);
            assign_block = completeBlock->getEmpty();
        }
        if (ret != un_comp->length) {
            memcpy(assign_block->data + assign_block->length, un_comp->data, ret * sizeof(char));
            assign_block->length += ret;
            completeBlock->backEmpty(assign_block);
            bam_number += find_divide_pos_and_get_read_number(assign_block).second;
            assign_block = completeBlock->getEmpty();

            memcpy(assign_block->data + assign_block->length, un_comp->data + ret,
                   (un_comp->length - ret) * sizeof(char));
            assign_block->length += (un_comp->length - ret);
            compress->backEmpty(un_comp);
            un_comp = compress->getUnCompressData();
            memcpy(assign_block->data + assign_block->length, un_comp->data, un_comp->length * sizeof(char));
            assign_block->length += un_comp->length;
            while (find_divide_pos(assign_block) != assign_block->length) {
                compress->backEmpty(un_comp);
                un_comp = compress->getUnCompressData();
                if (assign_block->length + un_comp->length > assign_block->data_size) {
                    change_data_size(assign_block);
                }
                memcpy(assign_block->data + assign_block->length, un_comp->data, un_comp->length * sizeof(char));
                assign_block->length += un_comp->length;
            }

            bam_number += find_divide_pos_and_get_read_number(assign_block).second;
            completeBlock->backEmpty(assign_block);
            assign_block = completeBlock->getEmpty();
        } else {
            if (ret != un_comp->length) {
                printf("ai nan ding\n");
                break;
            }
            assign_block->length += ret;
            last_use_block_length = 0;
            isclean = true;
            bam_number += un_comp->bam_number;
        }
    }
    if (assign_block->length != 0) {
        completeBlock->backEmpty(assign_block);
    }
    completeBlock->is_over();
    printf("Bam number is %d\n", bam_number);

}

void benchmark_bam_pack(BamCompleteBlock *completeBlock, samFile *output, sam_hdr_t *hdr) {

    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    hts_set_threads(output, 20);
    bam_complete_block *un_comp;
    long long ans = 0;
    long long res = 0;
    while (1) {
        un_comp = completeBlock->getCompleteBlock();
        if (un_comp == nullptr) {
            break;
        }
        int ret;
        while ((ret = (read_bam(un_comp, b, 0))) >= 0) {
            ans++;
        }
        res++;

        completeBlock->backEmpty(un_comp);
    }
    printf("Bam1_t Number is %lld\n", ans);
    printf("Block  Number is %lld\n", res);
}

void basic_status_pack(BamCompleteBlock *completeBlock, BamStatus *status) {
    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        printf("bam1_t is not ready\n");
    }
    bam_complete_block *un_comp;
    while (1) {
        un_comp = completeBlock->getCompleteBlock();
        if (un_comp == nullptr) {
            break;
        }
        int ret;
        while ((ret = (read_bam(un_comp, b, 0))) >= 0) {

            status->statusbam(b);
        }
        completeBlock->backEmpty(un_comp);
    }

}

int main(int argc, char *argv[]) {

    CLI::App app("RabbitBAM");

    string inputfile;

    string outputfile("./BAMStatus.html");
    int n_thread = 1;
    int n_thread_write = 1;
    int level = 6;

    bool is_tgs = false;


    CLI::App *benchmark_count = app.add_subcommand("benchmark_count", "read banchmark: count line");
    benchmark_count->add_option("-i", inputfile, "input File name")->required()->check(CLI::ExistingFile);
    benchmark_count->add_option("-o", outputfile, "output File name");
    benchmark_count->add_flag("--tgs", is_tgs, "Process TGS data");
    benchmark_count->add_option("-w,-@,-n,--threads", n_thread, "thread number");

    CLI::App *api_test = app.add_subcommand("api_test", "use api to read and write");
    api_test->add_option("-i", inputfile, "input File name")->required()->check(CLI::ExistingFile);
    api_test->add_option("-o", outputfile, "output File name");
    api_test->add_option("--nr", n_thread, "Read thread number");
    api_test->add_option("--nw", n_thread_write, "Write thread number, we suggest nr:nw = 1:4");
    api_test->add_flag("--tgs", is_tgs, "Process TGS data");
    api_test->add_option("-l,--level", level, "zip level");


    CLI11_PARSE(app, argc, argv);

    int total_threads = std::thread::hardware_concurrency();
    if (app.get_subcommands().size() > 1) {
        printf("you should input one command!!!\n");
        return 0;
    }
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "benchmark_count") == 0) {
        if(n_thread > total_threads - 4) n_thread = max(1, total_threads - 4);
        TDEF(fq)
        TSTART(fq)
        printf("Starting Running Benchmark Count\n");
        printf("BGZF_MAX_BLOCK_COMPLETE_SIZE is %d\n", BGZF_MAX_BLOCK_COMPLETE_SIZE);
        if (strcmp(outputfile.substr(outputfile.size() - 4).c_str(), "html") == 0) outputfile = ("./output.fastq");
        samFile *sin;
        sam_hdr_t *hdr;
        ofstream fout;
        fout.open(outputfile);
        if ((sin = sam_open(inputfile.c_str(), "r")) == NULL) {
            printf("Can`t open this file!\n");
            return 0;
        }
        if ((hdr = sam_hdr_read(sin)) == NULL) {
            return 0;
        }

        BamRead read(2000);
        BamCompress compress(8000, n_thread);
        BamCompleteBlock completeBlock(10);

        printf("Malloc Memory is Over\n");

        thread *read_thread = new thread(&read_pack_main, sin->fp.bgzf, &read);
        thread **compress_thread = new thread *[n_thread];

        for (int i = 0; i < n_thread; i++) {
            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(i, &cpuset);
            compress_thread[i] = new thread(&compress_pack_main, &read, &compress);
            int rc = pthread_setaffinity_np(compress_thread[i]->native_handle(), sizeof(cpu_set_t), &cpuset);

        }
        thread *assign_thread = new thread(&benchmark_pack, &compress, &completeBlock);
        read_thread->join();
        for (int i = 0; i < n_thread; i++) compress_thread[i]->join();
        assign_thread->join();
        sam_close(sin);
        printf("Wait num is %d\n", compress.wait_num);
        TEND(fq)
        TPRINT(fq, "time is : ");
    }
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "api_test") == 0) {
        if(n_thread + n_thread_write > total_threads - 6) {
            n_thread = max(1, (total_threads - 6) / 5);
            n_thread_write = max(1, total_threads - 6 - n_thread);
            printf("adjust read thread num to %d, write thread num to %d\n", n_thread, n_thread_write);
        }
        TDEF(fq)
        TSTART(fq)

        printf("Starting Running API Test\n");
        printf("BGZF_MAX_BLOCK_COMPLETE_SIZE is %d\n", BGZF_MAX_BLOCK_COMPLETE_SIZE);
        printf("output File Name is %s\n", outputfile.c_str());

//#define init_test
//#define debug_test 
//#define queue_test
#define parallel_test
//#define merge_parallel_test
//tagg


#ifdef init_test

        //int big_size = 800 << 10;
        int big_size = 256;
        //int big_size = 32 << 10;
        double t0 = GetTime();
        BamReader *reader = new BamReader(inputfile, n_thread);
        printf("new reader cost %lf\n", GetTime() - t0);

        t0 = GetTime();
        BamWriter *writer = new BamWriter(outputfile, reader->getHeader(), n_thread_write, level, big_size);
        printf("new writer cost %lf\n", GetTime() - t0);

        bam1_t *b;
        if ((b = bam_init1()) == NULL) {
            fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
        }
        long long num = 0;
        t0 = GetTime();
        double t1 = GetTime();
        while (reader->getBam1_t(b)) {
            num++;
            writer->write(b);
//            if (num%1000 == 0) printf("Bam1_t Num is %d\n",num);
        }
        printf("new process 1+2 cost %lf\n", GetTime() - t0);

        writer->over();
        printf("new total process cost %lf\n", GetTime() - t1);
#endif



#ifdef parallel_test


        if(is_tgs) {
            int big_size = 256;
            BamReader *reader = new BamReader(inputfile, n_thread, is_tgs);
            BamWriter *writer = new BamWriter(outputfile, reader->getHeader(), n_thread_write, level, big_size, is_tgs);
            bam1_t *b;
            if ((b = bam_init1()) == NULL) {
                fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
            }
            long long num = 0;
            while (reader->getBam1_t(b)) {
                num++;
                writer->write(b);
            }
            writer->over();
            cout << "Bam1_t Num : " << num << endl;
        } else {

            int big_size = 4 << 10;

            BamReader *reader = new BamReader(inputfile, n_thread, is_tgs);
            BamWriter *writer = new BamWriter(outputfile, reader->getHeader(), n_thread_write, level, big_size, is_tgs);

            const int vec_N = 40 << 10;
            std::vector < bam1_t * > b_vec[THREAD_NUM_P];
#pragma omp parallel for num_threads(THREAD_NUM_P) schedule(static)
            for (int i = 0; i < THREAD_NUM_P; i++) {
                for (int j = 0; j < vec_N; j++) {
                    bam1_t *item = bam_init1();
                    if (item == NULL) {
                        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
                    }
                    b_vec[i].push_back(item);
                }
            }

            long long num = 0;

            TDEF(inner)
            TSTART(inner)
            bool done = 0;
            while (done == 0) {
                auto res_vec = reader->getBam1_t_parallel(b_vec);
                int res_vec_size = res_vec.size();
                //printf("reader get vec size %d\n", res_vec_size);
                num += res_vec_size;
                if (res_vec_size == 0) break;
                writer->write_parallel(res_vec);
            }
            writer->over_parallel();
            TEND(inner)
            TPRINT(inner, "inner time is : ")
            cout << "Bam1_t Num : " << num << endl;
        }
#endif
        TEND(fq)
        TPRINT(fq, "time is : ");


    }
}



