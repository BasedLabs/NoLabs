# Base image
FROM python:3.11-buster

ENV POETRY_VERSION=1.6.1
ENV POETRY_VIRTUALENVS_CREATE=false
ENV PYTHONUNBUFFERED=1

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    && rm -rf /var/lib/apt/lists/*

RUN pip install poetry
ENV PYTHONBUFFERED=1

WORKDIR /app
COPY . /app

RUN poetry install --no-root --directory /app/nolabs

ENTRYPOINT ["poetry", "run", "--directory", "/app/nolabs"]