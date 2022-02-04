#!/bin/bash

inputs=./tests/*.json

frames=()
outputext=png

exebase=bomat
outdir=./tests/
expectedoutdir=./tests/expected-output
use_stdin="false"
use_pushpop="false"

source ./submodules/bat/test.sh

