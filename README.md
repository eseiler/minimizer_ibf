# Minimizer example and benchmark
This repository contains reference data and (minimal) code for computing minimizers and (soon) the IBF.

## Requirements
* cmake >= 3.2
* g++ >= 7

## Setup
```console
git clone --recursive https://github.com/eseiler/minimizer_ibf
cd minimizer_ibf
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make minimizer_example
make minimizer_benchmark
./minimizer_example
./minimizer_benchmark
```
Note the `--recursive` when cloning.

If the correct compiler is not detected by default, pass `-DCMAKE_CXX_COMPILER=executable` to the `cmake` command. `executable` refers to the compiler, e.g. `g++-9`.

See `./minimizer_benchmark --help` for more options on how to run benchmarks.
