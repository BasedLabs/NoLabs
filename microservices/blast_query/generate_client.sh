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

# Check if the EMAIL argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 EMAIL"
    exit 1
fi

EMAIL=$1

echo 'Installation: https://openapi-generator.tech/docs/installation'

DOCKER_IMAGE_NAME="blast_query"

# Build the Docker image
docker build -t $DOCKER_IMAGE_NAME -f build/Dockerfile .

# Run the Docker container in the background, passing in the EMAIL
docker run -d --name blast_query -p 5743:5743 -e EMAIL=$EMAIL $DOCKER_IMAGE_NAME --host=0.0.0.0 --port=5743

# Generate the Python client using OpenAPI Generator
npx @openapitools/openapi-generator-cli generate \
    -i http://127.0.0.1:5743/openapi.json \
    -g python \
    -o ./client \
    --additional-properties=packageName=blast_query_microservice

echo 'Use pip install ./client'