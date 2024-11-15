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

FROM ghcr.io/basedlabs/rfdiffusion-base:1.0.0

ENV RFDIFFUSION_PATH=/app/RFdiffusion

ENV HYDRA_FULL_ERROR=1

WORKDIR /app
RUN mkdir $RFDIFFUSION_PATH/models
COPY service /app
COPY --from=requirements /app/requirements.txt /app/requirements.txt
RUN pip install -r requirements.txt

ENTRYPOINT []