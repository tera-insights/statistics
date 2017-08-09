#!/bin/bash -e 

mkdir -p test/bin
clang++ -std=c++11 test/memory-conscious-hashing.cc -o test/bin/test-mch