# RabbitBAM

An ultra-fast I/O framework for BAM files.

## Dependancy

- gcc 8.5.0 or newer

## Installation

```bash
git clone https://github.com/RabbitBio/RabbitBAM
cd RabbitBAM
bash build_all.sh
export LD_LIBRARY_PATH=./ext/libdeflate-1.20/build:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=./ext/htslib-1.20:$LD_LIBRARY_PATH
```

## Benchmark

refer to `./rabbitbam -h`

