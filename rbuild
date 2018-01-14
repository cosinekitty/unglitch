#!/bin/bash
mkdir -p release
pushd release > /dev/null
cmake -DCMAKE_BUILD_TYPE=Release .. && make && make install
popd
cp -v bin/Release/unglitch ~/bin/unglitch
