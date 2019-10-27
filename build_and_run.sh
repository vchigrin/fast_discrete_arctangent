#!/bin/sh
clang++ -std=c++11 -O2 -g -Wl,-lgtest_main -Wl,-lgtest -Wl,-lpthread fast_discrete_arctangent_test.cc && ./a.out
