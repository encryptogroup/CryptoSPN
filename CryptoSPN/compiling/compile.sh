#!/bin/sh

original_dir=$PWD
tmp_dir=$(mktemp -d -t cryptospn-XXXXXXXXXX --tmpdir="${2}/src/examples/")
cp "${3}" "${tmp_dir}"/
cp -r "${1}"/aby_files/selection_blocks "${tmp_dir}"/selection_blocks
cp "${5}" "${tmp_dir}"/CMakeLists.txt
mkdir "${tmp_dir}"/build
cd "${tmp_dir}"/build || exit 1
cmake ..  -DCMAKE_BUILD_TYPE=Release
# workaround for relic issues when you build your own application (requires ABY to have been built before):
make >/dev/null 2>&1
cp "${2}"/build/lib/librelic_s.a build/lib/ || exit 2
make
status="$?"
sleep "1s"
cp "${4}" "${original_dir}"
rm -rf "${tmp_dir}"
exit $status