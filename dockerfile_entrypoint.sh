#!/bin/bash

npm --prefix /app/src/server/frontend run dev -- --host &
python src/server/run.py --host=127.0.0.1 $TEST