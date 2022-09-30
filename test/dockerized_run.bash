#!/bin/bash

set -x
set -e

mkdir -p results
mkdir -p data

DOCKERFILE_TAG="dockerized_smk"

# Make sure snakemake.smk is in external/ and requirements.txt is updated
# cp ../snakefile.smk external/
cp ../pyproject.toml .

docker build . --tag ${DOCKERFILE_TAG} --file DOCKERFILE
# --no-cache 

# --rm removes the container after it exits
# -v mounts a volume into the container
# --it makes the container interactive
time docker run -it --rm -v ${PWD}/external:/data/external -v ${PWD}/mzml:/data/mzml -v $PWD/bin:/data/bin -v $PWD/results:/data/results ${DOCKERFILE_TAG}