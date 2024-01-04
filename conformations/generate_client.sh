echo 'Installation: https://openapi-generator.tech/docs/installation'
sudo docker run --name conformations --gpus all --net=host conformations --host=0.0.0.0 --port=5731 &
npx @openapitools/openapi-generator-cli generate -i http://127.0.0.1:5731/openapi.json -g python -o ./client --additional-properties=packageName=conformations_microservice
echo 'Use pip install conformations/client'