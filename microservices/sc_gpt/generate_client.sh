#!/bin/bash

#
# IMPORTANT: run this file like this:
# $ ./generate_client.sh EMAIL
#

# Check for Docker installation
if ! command -v docker &> /dev/null; then
    echo "Docker could not be found. Please install Docker."
    exit 1
fi

# Check for npx installation
if ! command -v npx &> /dev/null; then
    echo "npx could not be found. Please install Node.js and npm."
    exit 1
fi


echo 'Installation: https://openapi-generator.tech/docs/installation'

DOCKER_IMAGE_NAME="sc_gpt"

# Build the Docker image
docker build -t $DOCKER_IMAGE_NAME -f build/Dockerfile .

# Run the Docker container in the background
docker run -d --name sc_gpt -p 5744:5744 sc_gpt --host=0.0.0.0 --port=5744

# Generate the Python client using OpenAPI Generator
npx @openapitools/openapi-generator-cli generate \
    -i http://127.0.0.1:5744/openapi.json \
    -g python \
    -o ./client \
    --additional-properties=packageName=sc_gpt_microservice

echo 'Use pip install ./client'