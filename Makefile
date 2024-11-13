# Load environment variables from .env file
ifneq (,$(wildcard ./.env))
    include .env
    export
endif

flower:
	celery --broker=${CELERY_BROKER_URL} flower --port=5555
install-openapi-generator:
	npm install -g openapi-typescript-codegen
generate-client:
	@openapi --input 'http://127.0.0.1:8000/openapi.json' --output frontend/src/api/client --client axios
install-mock-server:
	npm install -g @stoplight/prism-cli
start-mock-server:
	prism mock http://127.0.0.1:${UVICORN_HOST}/openapi.json
prepare-nolabs:
	scripts/prepare-nolabs.sh
download-diffdock-weights:
	@echo "Downloading model weights..."
	python3 microservices/diffdock/scripts/download_weights.py
	@echo "Download complete!"
download-esmfold-weights:
	@echo "Downloading model weights..."
	pip3 install transformers[torch] --verbose
	python3 microservices/esmfold/scripts/download_weights.py
	@echo "Download complete!"