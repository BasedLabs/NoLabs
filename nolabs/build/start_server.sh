npm --prefix frontend run dev -- --host &
NOLABS_ENVIRONMENT=dev poetry run uvicorn nolabs.application.api:app --host=0.0.0.0 --port=8000 --workers 4