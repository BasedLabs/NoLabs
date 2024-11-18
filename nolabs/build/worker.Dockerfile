# Generate workable requirements.txt from Poetry dependencies
FROM python:3.11-slim as requirements

RUN pip install poetry poetry-plugin-export

ENV POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1 \
    POETRY_CACHE_DIR=/tmp/poetry_cache

WORKDIR /app
COPY nolabs/pyproject.toml /app
RUN poetry lock --no-update
RUN poetry export -f requirements.txt --without-hashes -o requirements.txt

FROM python:3.11-buster
COPY --from=requirements /app/requirements.txt /app/requirements.txt
RUN pip install -r /app/requirements.txt
COPY nolabs /app/nolabs
WORKDIR /app/nolabs
