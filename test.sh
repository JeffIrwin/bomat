#!/bin/bash

#=======================================================================

#inputs=./tests/5.json
inputs=./tests/*.json

frames=()
outputext=png

exebase=bomat
#exebase=$(grep -P -o 'set\(\s*PROJECT\s*\K.*(?=\s*\))' CMakeLists.txt)

outdir=./tests/
expectedoutdir=./tests/expected-output
use_stdin="false"
use_pushpop="false"

# The parentheses are necessary to stop the tests from exiting on success?
# Weird test but ok
(
	source ./submodules/bat/test.sh
)
if [[ "$?" != "0" ]]; then
	exit -1
fi

#=======================================================================

# Dry run syntax check on all examples.  Can this be built into bat?

examples=./examples/*.json
for e in ${examples}; do
	#echo ${e}

	# TODO: handle errors
	./build/bomat -d "${e}" || exit -2

done

#=======================================================================

echo "All tests and lints passed!"
echo ""

#=======================================================================

