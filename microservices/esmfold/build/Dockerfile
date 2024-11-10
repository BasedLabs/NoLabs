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

# Create a Conda environment with Python 3.11
RUN conda create -n myenv python=3.10 && conda clean -a -y

# Activate the environment and install packages
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]
RUN pip install poetry poetry-plugin-export

# Set up poetry and dependencies
WORKDIR /app
COPY pyproject.toml /app
COPY --from=requirements /app/requirements.txt /app/requirements.txt
RUN pip install -r requirements.txt

# Copy application files
COPY service /app/service
WORKDIR /app/service

# Run model loading step within the environment
RUN python -c "from transformers import EsmForProteinFolding; EsmForProteinFolding.from_pretrained('facebook/esmfold_v1', low_cpu_mem_usage=True)"
