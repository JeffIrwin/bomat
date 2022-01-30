#!/bin/bash

inputs=./examples/*.json

frames=()
outputext=png

exebase=bomat
outdir=./examples/
expectedoutdir=./examples/expected-output
use_stdin="false"
use_pushpop="false"

# HACK
source ./submodules/bat/test.sh || exit 0

