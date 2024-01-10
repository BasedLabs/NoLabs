sudo docker build --progress=plain -t base_conformations -f build/base.Dockerfile .
DOCKER_IMAGE=$(docker images -q base_conformations)
docker tag $DOCKER_IMAGE ghcr.io/basedlabs/conformations-base:latest
docker push ghcr.io/basedlabs/conformations-base:latest