#!/bin/sh

set -e
type cmake3 >/dev/null 2>&1 && alias cmake=cmake3

# Every path is quoted: on macOS the source tree may sit under
# "~/Library/Mobile Documents/..." (iCloud Drive), and an unquoted
# expansion splits that path at the space.
main_dir=$(dirname "$(readlink -f "$0")")
src_dir=$main_dir/src
linkdef_dir=$main_dir/linkdef
build_dir=$main_dir/.build
bin_dir=$main_dir/bin

##### macOS
if [ "$(uname)" = 'Darwin' ] && [ ! -e "$src_dir/Dict_rdict.pcm" ]; then
    rootcling -f "$src_dir/Dict.cc" -c TVector3.h TParticle.h \
	      "$linkdef_dir/LinkDef.h"
fi

mkdir -p "$build_dir"
cd "$build_dir"
cmake .. -DCMAKE_INSTALL_PREFIX="$main_dir" \
      -DCMAKE_INSTALL_RPATH_USE_LINK_PATH="ON"
cmake --build .
cmake --install .

if [ "$(uname)" = 'Darwin' ] && [ ! -e "$bin_dir/Dict_rdict.pcm" ]; then
    cp "$src_dir/Dict_rdict.pcm" "$bin_dir"
fi
