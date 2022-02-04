#!/bin/bash

pushd .

date

time ./build/bomat || exit -1

popd

