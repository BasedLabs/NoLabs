npm --prefix frontend run dev -- --host &
NOLABS_ENVIRONMENT=dev poetry run uvicorn nolabs.api:app --host=0.0.0.0 --port=8000