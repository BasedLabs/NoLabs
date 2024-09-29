flower:
	celery --broker=redis://localhost:6379/0 flower --port=5555
install-openapi-generator:
	npm install -g openapi-typescript-codegen
generate-client:
	@openapi --input 'http://127.0.0.1:8000/openapi.json' --output frontend/src/api/client --client axios