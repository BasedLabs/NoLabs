if [ ! -d "$RFdiffusion" ]; then
  git clone https://github.com/RosettaCommons/RFdiffusion.git
fi
cd RFdiffusion
docker buildx build --progress=plain -t rfdiffusion -f docker/Dockerfile .

DOCKER_IMAGE=$(docker images -q rfdiffusion)
docker tag $DOCKER_IMAGE ghcr.io/basedlabs/rfdiffusion-base:1.0.0
docker push ghcr.io/basedlabs/rfdiffusion-base:1.0.0