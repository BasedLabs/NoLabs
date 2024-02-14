# Use a base image with CUDA support that matches the required PyTorch and CUDA versions
FROM pytorch/pytorch:1.11.0-cuda11.3-cudnn8-runtime

# Install necessary tools and libraries
RUN apt-get update && apt-get install -y wget git

# Set a working directory inside the container
WORKDIR /app

# Copy the requirements file and install Python dependencies
COPY ./requirements.txt /app/
# Ensure pip is up-to-date
RUN pip install --upgrade pip
# Install dependencies from requirements.txt
RUN pip install --default-timeout=100 -r requirements.txt

RUN apt-get update && apt-get install -y build-essential

# Install PyTorch Geometric and its dependencies
RUN pip install torch==1.11.0 torchaudio torchvideo
# It's important to install these after PyTorch is installed and to use the compatible CUDA version
RUN pip install torch-geometric==2.0.4
RUN pip install --no-index torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-1.11.0+cu113.html

# Install additional dependencies for DiffDock
RUN pip install PyYAML scipy "networkx[default]" biopython rdkit-pypi e3nn spyrmsd pandas biopandas
RUN pip install "fair-esm[esmfold]"
RUN pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'

# Clone the DiffDock repository into the specified directory structure
RUN mkdir -p diffdock/diffdock/diffdock_src && git clone https://github.com/gcorso/DiffDock.git diffdock/diffdock/diffdock_src

# Copy the application files into the container
COPY ./diffdock /app/diffdock

ENTRYPOINT uvicorn diffdock.api:app --host $HOST --port $PORT