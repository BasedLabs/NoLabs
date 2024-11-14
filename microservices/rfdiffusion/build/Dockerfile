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

FROM ghcr.io/basedlabs/rfdiffusion:latest

ENV HYDRA_FULL_ERROR=1
ARG rfdiffusion_path=/app/protein_design/RFdiffusion

WORKDIR /app/protein_design
WORKDIR /app
RUN apt-get update
RUN apt-get -y install wget git
RUN git clone https://github.com/RosettaCommons/RFdiffusion $rfdiffusion_path
RUN mkdir $rfdiffusion_path/models
RUN ls -l $rfdiffusion_path
WORKDIR $rfdiffusion_path/scripts
RUN chmod +x download_models.sh
RUN sed -i 's/\r$//' download_models.sh # Fix line endings to UNIX
RUN ./download_models.sh $rfdiffusion_path/models
WORKDIR /app
ADD . /app
ADD requirements.txt /app
RUN pip install -r requirements.txt
RUN ls -l

ENTRYPOINT ["uvicorn", "protein_design.api:app", "--workers", "2"]