# Minimizer IBF
This repository contains an implementation of the minimizer IBF.

## Requirements
* cmake >= 3.2
* g++ >= 7

## Setup
```console
git clone https://github.com/eseiler/minimizer_ibf
cd minimizer_ibf
git checkout seqan3
git submodule update --init --force --recursive
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-march=native
make build_min_ibf
make search_min_ibf
```

```console
# Example data
cp ../data/ibf.tar.lzma .
tar xf ibf.tar.lzma
# Build, use 1 MiB per bin, otherwise use defaults
./build_min_ibf ibf/bins/ ibf.out --bits 1048576
# Search for all reads that originate from bin 0
./search_min_ibf ibf/reads/bin_00.fastq ibf.out
./search_min_ibf ibf/reads/bin_00.fastq ibf.out --error 1
./search_min_ibf ibf/reads/bin_00.fastq ibf.out --error 2
```

If the correct compiler is not detected by default, pass `-DCMAKE_CXX_COMPILER=executable` to the `cmake` command. `executable` refers to the compiler, e.g. `g++-9`.

See `./build_min_ibf --help` or `./search_min_ibf --help` for more options on how to run benchmarks.
