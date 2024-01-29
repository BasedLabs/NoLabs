npm --prefix frontend run dev -- --host &
echo $1
uvicorn nolabs.api:app --host=0.0.0.0 $1