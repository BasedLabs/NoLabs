# Generate workable requirements.txt from Poetry dependencies
FROM python:3.10-slim as requirements

RUN pip install poetry poetry-plugin-export

ENV POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1 \
    POETRY_CACHE_DIR=/tmp/poetry_cache

WORKDIR /app
COPY pyproject.toml /app
RUN poetry lock --no-update
RUN poetry export -f requirements.txt --without-hashes -o requirements.txt

FROM pytorch/pytorch:1.13.1-cuda11.6-cudnn8-runtime
RUN apt-get update && apt-get install -y wget git build-essential
WORKDIR /app
COPY --from=requirements /app/requirements.txt /app/requirements.txt
RUN pip install -r /app/requirements.txt --verbose
RUN git clone https://github.com/dauparas/ProteinMPNN.git /app/ProtMPNN
COPY service /app

ENTRYPOINT []