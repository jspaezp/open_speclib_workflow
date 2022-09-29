#!/bin/bash

# This takes a while so expect to wait ~30 minutes for it to compile


if [ ! -d pwiz ] ;  then
   git clone --depth 1 https://github.com/ProteoWizard/pwiz.git
else
   cd pwiz
   git pull
   cd ..
fi

BUILD_COMMAND="bash quickbuild.sh -j4 --abbreviate-paths --hash optimization=space address-model=64 pwiz_tools/BiblioSpec toolset=gcc"

if [ -z "$(which docker)" ] ; then
   echo "Docker not found. Building locally in 5 seconds."
   sleep 5
   bash -c "cd pwiz && ${BUILD_COMMAND}"
else
   docker build . --tag my-gcc --file DOCKERFILE
   time docker run -it --rm -v $PWD:/build -w /build/pwiz my-gcc /bin/bash -c "${BUILD_COMMAND}"
fi

cp pwiz/build-linux-x86_64/BiblioSpec/bibliospec-bin-linux-*-release-*.tar.bz2 .
tar -x --bzip2 -f bibliospec-bin-linux-*-release-*.tar.bz2
mkdir -p ../../test/bin
cp Blib* ../../test/bin/.