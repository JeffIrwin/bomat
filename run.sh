#!/bin/bash

pushd .

date

time ./build/bomat.exe || exit -1

popd

