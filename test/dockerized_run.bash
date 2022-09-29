
mkdir dockervol

DOCKERFILE_TAG="dockerized_smk"

docker build . --tag ${DOCKERFILE_TAG} --file DOCKERFILE

# --rm removes the container after it exits
# -v mounts a volume into the container
# --it makes the container interactive
time docker run -it --rm -v $PWD/bin:/data/bin -v $PWD/dockervol:/data/results ${DOCKERFILE_TAG}