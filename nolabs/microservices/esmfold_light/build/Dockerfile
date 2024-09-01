FROM python:3.11-alpine

RUN pip install poetry

ENV POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1 \
    POETRY_CACHE_DIR=/tmp/poetry_cache

WORKDIR /app

COPY pyproject.toml poetry.lock ./

RUN poetry install && rm -rf $POETRY_CACHE_DIR

COPY service /app/service

WORKDIR /app/service

ENTRYPOINT ["poetry", "run", "python", "api.py"]