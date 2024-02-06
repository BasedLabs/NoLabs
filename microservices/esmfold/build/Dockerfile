# Use a PyTorch base image with CUDA support
FROM pytorch/pytorch:1.13.1-cuda11.6-cudnn8-runtime

# Install Python 3.8 and curl
RUN apt-get update && apt-get install -y python3.8 curl python3-pip python3.8-venv

# Set a working directory
WORKDIR /app

# Copy only the requirements.txt first to leverage Docker cache
COPY ./requirements.txt /app/

# Install dependencies
RUN python3.8 -m pip install --upgrade pip
RUN python3.8 -m pip install --default-timeout=100 -r requirements.txt

# Pre-download ESM-Fold model weights
RUN python3.8 -c "from transformers import EsmForProteinFolding; EsmForProteinFolding.from_pretrained('facebook/esmfold_v1', low_cpu_mem_usage=True)"

# Now copy the rest of the application files into the container
COPY ./esmfold /app/esmfold

# Set the command to run the FastAPI app
ENTRYPOINT ["uvicorn", "esmfold.api:app"]
