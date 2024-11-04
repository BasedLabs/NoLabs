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

RUN apt-get update; apt-get install curl gpg -y; \
mkdir -p /etc/apt/keyrings; \
curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg; \
echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_20.x nodistro main" | tee /etc/apt/sources.list.d/nodesource.list; \
 apt-get update && apt-get install -y nodejs;

RUN pip install poetry

COPY . /app
WORKDIR /app/frontend
RUN npm install
RUN npm run build

ENV PYTHONBUFFERED=1
ENV COMPOSE=1

RUN poetry install --no-root
WORKDIR /app
USER root
RUN chmod +x build/start_server.sh

RUN pip install uvicorn

RUN ls -la
ENTRYPOINT ["sh", "./build/start_server.sh"]