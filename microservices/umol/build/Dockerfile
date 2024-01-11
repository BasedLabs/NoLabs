FROM nvcr.io/nvidia/jax:23.08-t5x-py3

# Install Python 3.8
RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential \
  git \
  curl \
  vim \
  python3-dev \
  python3-pip

# Download the model weights and place them in the specified directory
RUN mkdir -p /app/umol/umol_source
RUN curl -o /app/umol/umol_source/params.pkl -L "https://huggingface.co/thomasshelby/uMol_params/resolve/main/params.pkl?download=true"

# Set a working directory
WORKDIR /app/umol
WORKDIR /app

# Copy the application files into the container
COPY ./requirements.txt /app

# Install dependencies
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

RUN pip install --upgrade "jax[cuda11_pip]==0.4.13" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

COPY ./umol /app/umol

RUN ls -l /app/umol
RUN ls -l /app/umol/umol_source

# Set the command to run the FastAPI app
ENTRYPOINT uvicorn umol.api:app --host $HOST --port $PORT
