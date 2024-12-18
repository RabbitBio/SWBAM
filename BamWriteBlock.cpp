#include "BamWriteBlock.h"

BamWriteBlockConfig::BamWriteBlockConfig() {}

BamWriteBlockConfig::BamWriteBlockConfig(int Buffer_number) {
    this->Buffer_number = Buffer_number;
    this->write_number = Buffer_number + 10;
    this->complete = 0;
}

BamWriteBlock::BamWriteBlock() {};

BamWriteBlock::BamWriteBlock(BamWriteBlockConfig *config){
    this->config = config;
    this->buffer = new bam_block*[this->config->Buffer_number];
    for (int i = 0; i<this->config->Buffer_number; ++i) this->buffer[i] = new bam_block;
    this->compress = new int[this->config->write_number];
    this->compress_bg = 0;
    this->compress_ed = 0;
    this->read = new int[this->config->write_number];
    this->read_bg = 0;
    this->read_ed = this->config->Buffer_number;
    for (int i = this->read_bg; i < this->read_ed; i++) this->read[i] = i;
}

std::pair<bam_block *, int> BamWriteBlock::getEmpty() {
    while (read_bg == read_ed) {
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        //this_thread::yield();
    }
    int num = read_bg;
    read_bg = (read_bg + 1) % config->write_number;
    return std::pair<bam_block *, int>(this->buffer[read[num]], read[num]);

}

void BamWriteBlock::inputblock(int id) {
    compress[compress_ed] = id;
    compress_ed = (compress_ed + 1) % config->write_number;
}

std::pair<bam_block *, int> BamWriteBlock::getCompressdata() {
    mtx_compress.lock();
    while (compress_ed == compress_bg) {
        mtx_compress.unlock();
        std::this_thread::sleep_for(std::chrono::milliseconds(3));
        //this_thread::yield();
        if (config->complete) return std::pair<bam_block *, int>(NULL, -1);
        mtx_compress.lock();
    }
    int num = compress_bg;
    compress_bg = (compress_bg + 1) % config->write_number;
    mtx_compress.unlock();
    return std::pair<bam_block *, int>(this->buffer[compress[num]], compress[num]);
}

void BamWriteBlock::backempty(int id) {
    mtx_read.lock();
    read[read_ed] = id;
    read_ed = (read_ed + 1) % config->write_number;
    mtx_read.unlock();
}

bool BamWriteBlock::isComplete() {
    return config->complete;
};

void BamWriteBlock::ReadComplete() {
    config->complete = 1;
}
