
echo 'Login here https://stackoverflow.com/questions/75924590/pushing-to-github-container-registry-unauthorized-unauthenticated-user-cannot'

DOCKER_IMAGE=$(docker images -q rfdiffusion)
docker tag $DOCKER_IMAGE ghcr.io/basedlabs/rfdiffusion:latest
docker push ghcr.io/basedlabs/rfdiffusion:latest

#if [ ! -d "$RFdiffusion" ]; then
#  git clone https://github.com/RosettaCommons/RFdiffusion.git
#fi
#cd RFdiffusion
#docker buildx build --progress=plain -t rfdiffusion -f docker/Dockerfile .