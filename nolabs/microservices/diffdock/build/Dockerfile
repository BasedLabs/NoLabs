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

FROM pytorch/pytorch:1.11.0-cuda11.3-cudnn8-runtime

RUN apt-get update && apt-get install -y wget git build-essential

WORKDIR /app

# Clone the DiffDock repository and checkout the specified commit
RUN git clone https://github.com/gcorso/DiffDock.git && \
    cd DiffDock && \
    git checkout a6c5275

# Clone the ESM repository, check out the specified commit, and install
RUN if [ ! -d "DiffDock/esm" ]; then \
        cd DiffDock && \
        git clone https://github.com/facebookresearch/esm && \
        cd esm && \
        git checkout ca8a710 && \
        pip install -e .; \
    fi

RUN wget https://sourceforge.net/projects/smina/files/smina.static/download -O /app/DiffDock/smina && \
    chmod +x /app/DiffDock/smina

# Download and install gnina
RUN wget https://github.com/gnina/gnina/releases/download/v1.0.3/gnina -O /app/DiffDock/gnina && \
    chmod +x /app/DiffDock/gnina

COPY --from=requirements /app/requirements.txt /app/requirements.txt

RUN pip install -r /app/requirements.txt

ADD service /app/service

WORKDIR /app/service
