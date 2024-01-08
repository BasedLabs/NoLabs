# CUDA 11 base image
ARG CUDA="11.4.3"
ARG CUDNN="8"
ARG TAG="devel"
ARG OS="ubuntu20.04"
FROM nvidia/cuda:${CUDA}-cudnn${CUDNN}-${TAG}-${OS}

# Install Python 3.8
RUN apt-get update && apt-get install -y python3.8 python3-pip curl

# Set a working directory
WORKDIR /app

# Copy the application files into the container
COPY ./esmfold /app/esmfold
COPY ./requirements.txt /app

# Install dependencies
RUN python3.8 -m pip install --upgrade pip
RUN python3.8 -m pip install -r requirements.txt

# Pre-download ESM-Fold model weights
RUN python3.8 -c "from transformers import EsmForProteinFolding; EsmForProteinFolding.from_pretrained('facebook/esmfold_v1', low_cpu_mem_usage=True)"

# Set the command to run the FastAPI app
ENTRYPOINT uvicorn esmfold.api:app --host $HOST --port $PORT
