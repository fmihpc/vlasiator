#!/bin/bash

set -e   # Abort on error

# Header-only / source dependencies that used to be git submodules.
# Cloned into ./submodules/, matching the paths expected by the Makefile.
# Pinned to the same commits the git submodules last pointed at, so
# updating these versions is a deliberate, reviewable change.

FSGRID_COMMIT="34f8ffbe72db0fe3119d9013a21bc4df117b5307"
DCCRG_COMMIT="f086044ff1cf683ca125e61c9e43e4584b6a00bf"
EIGEN_COMMIT="3147391d946bb4b6c68edd901f2add6ac1f31f8c"
VECTORCLASS_COMMIT="f4617df57e17efcd754f5bbe0ec87883e0ed9ce6"
VECTORCLASS_ADDON_COMMIT="600d67becf8144cff4da77b47dd832ce9a98ab2e"
HASHINATOR_COMMIT="4ec66208bb5d820e766cfa0765a903ba5de4e4ca"

git_use_commit() {
	if [[ z$1 != "z" ]]; then
		git fetch origin "$1"
		git checkout "$1"
	fi
}

mkdir -p submodules
cd submodules

rm -rf fsgrid
git clone https://github.com/fmihpc/fsgrid.git
cd fsgrid
git_use_commit "$FSGRID_COMMIT"
cd ..

rm -rf dccrg
git clone -b vlasiator-version https://github.com/fmihpc/dccrg.git
cd dccrg
git_use_commit "$DCCRG_COMMIT"
cd ..

rm -rf eigen
git clone -b master https://gitlab.com/libeigen/eigen.git
cd eigen
git_use_commit "$EIGEN_COMMIT"
cd ..

rm -rf vectorclass
git clone https://github.com/vectorclass/version2 vectorclass
cd vectorclass
git_use_commit "$VECTORCLASS_COMMIT"
cd ..

rm -rf vectorclass-addon
git clone https://github.com/vectorclass/add-on vectorclass-addon
cd vectorclass-addon
git_use_commit "$VECTORCLASS_ADDON_COMMIT"
cd ..

rm -rf hashinator
git clone https://github.com/fmihpc/hashinator.git
cd hashinator
git_use_commit "$HASHINATOR_COMMIT"
cd ..
