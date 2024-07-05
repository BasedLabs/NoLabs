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

DOCKER_IMAGE_NAME="external_data_query"

sudo docker build -t $DOCKER_IMAGE_NAME -f build/Dockerfile .
# Run the Docker container in the background
# Replace the Docker image name with the appropriate one
sudo docker run -d --name external_data_query -p 5739:5739 external_data_query --host=0.0.0.0 --port=5739

# Generate the Python client using OpenAPI Generator
npx @openapitools/openapi-generator-cli generate \
    -i http://127.0.0.1:5739/openapi.json \
    -g python \
    -o ./client \
    --additional-properties=packageName=external_data_query_microservice

echo 'Use pip install ./client'
