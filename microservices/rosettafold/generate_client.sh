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

# Generate the Python client using OpenAPI Generator
npx @openapitools/openapi-generator-cli generate \
    -i http://127.0.0.1:5738/openapi.json \
    -g python \
    -o ./client \
    --additional-properties=packageName=rosettafold_microservice

echo 'Use pip install ./client'
