#!/bin/bash

#
# IMPORTANT: run this file like this:
# $ ./generate_client.sh YOUR_OPENAI_API_KEY_HERE
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

# Check if the OPENAI_API_KEY argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 OPENAI_API_KEY"
    exit 1
fi

OPENAI_API_KEY=$1

echo 'Installation: https://openapi-generator.tech/docs/installation'

DOCKER_IMAGE_NAME="biobuddy"

# Build the Docker image
sudo docker build -t $DOCKER_IMAGE_NAME -f build/Dockerfile .

# Run the Docker container in the background, passing in the OPENAI_API_KEY
sudo docker run -d --name biobuddy -p 5738:5738 -e OPENAI_API_KEY=$OPENAI_API_KEY $DOCKER_IMAGE_NAME --host=0.0.0.0 --port=5738

# Generate the Python client using OpenAPI Generator
npx @openapitools/openapi-generator-cli generate \
    -i http://127.0.0.1:5738/openapi.json \
    -g python \
    -o ./client \
    --additional-properties=packageName=biobuddy_microservice

echo 'Use pip install ./client'