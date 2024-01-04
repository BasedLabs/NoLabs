docker run --rm \
    -v $PWD:/local openapitools/openapi-generator-cli generate \
    -i http://127.0.0.1:5731/openapi.json \
    -g python \
    -o python-client.py
