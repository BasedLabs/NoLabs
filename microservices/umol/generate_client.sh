#!/bin/bash

# Make sure Docker is installed
if ! command -v docker &> /dev/null
then
    echo "Docker could not be found. Please install Docker."
    exit
fi

# Make sure npx is installed
if ! command -v npx &> /dev/null
then
    echo "npx could not be found. Please install Node.js and npm."
    exit
fi

echo 'Installation: https://openapi-generator.tech/docs/installation'

# Specify the Docker image to use (replace with umol_app_cuda11 or umol_app_cuda12 as needed)
DOCKER_IMAGE_NAME="umol" # or "umol_app_cuda12"

sudo docker build -t $DOCKER_IMAGE_NAME -f build/Dockerfile .

# Run the Docker container in the background
# Replace the Docker image name with the appropriate one
sudo docker run -d --name umol --gpus all -e HOST=0.0.0.0 -e PORT=5735 -p 5735:5735 $DOCKER_IMAGE_NAME

# Generate the Python client using OpenAPI Generator
npx @openapitools/openapi-generator-cli generate \
    -i http://127.0.0.1:5735/openapi.json \
    -g python \
    -o ./client \
    --additional-properties=packageName=umol_microservice

echo 'Use pip install ./client'

# Note: The Docker container named 'umol' will continue to run. Stop and remove it if necessary.
