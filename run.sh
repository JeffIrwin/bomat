#!/bin/bash

date

# 1st pass: calculate eigenvalues
time ./build/bomat tests/1.json -e

# 2nd pass: load eigenvalues from 1st pass and draw plot
time ./build/bomat tests/1.json -p

