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

# Run the Docker container in the background
# Replace 'p2rank_app' with the actual name of your Docker image
# Adjust the port to match the port you're using for the FastAPI application
sudo docker run -d --name p2rank -e HOST=0.0.0.0 -e PORT=5731 -p 5731:5731 p2rank_app

# Generate the Python client using OpenAPI Generator
npx @openapitools/openapi-generator-cli generate \
    -i http://127.0.0.1:5731/openapi.json \
    -g python \
    -o ./client \
    --additional-properties=packageName=p2rank_microservice

echo 'Use pip install ./client'
