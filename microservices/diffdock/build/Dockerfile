# Use a base image with CUDA support that matches the required PyTorch and CUDA versions
FROM pytorch/pytorch:1.11.0-cuda11.3-cudnn8-runtime

# Install necessary tools and libraries
RUN apt-get update && apt-get install -y wget git

# Set up Anaconda and create a conda environment for DiffDock
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh
RUN bash /miniconda.sh -b -p /miniconda
ENV PATH="/miniconda/bin:${PATH}"
RUN conda init bash
RUN . ~/.bashrc
RUN conda create --name diffdock python=3.9 -y
SHELL ["conda", "run", "-n", "diffdock", "/bin/bash", "-c"]

# Install PyTorch, CUDA, and PyTorch Geometric
RUN conda install pytorch==1.11.0 pytorch-cuda=11.3 -c pytorch -c nvidia -y
RUN pip install torch-scatter torch-sparse torch-cluster torch-spline-conv torch-geometric==2.0.4 -f https://data.pyg.org/whl/torch-1.11.0+cu113.html

# Set a working directory inside the container
WORKDIR /app

# Clone the DiffDock repository into the specified directory structure
RUN mkdir -p diffdock/diffdock/diffdock_src && git clone https://github.com/gcorso/DiffDock.git diffdock/diffdock/diffdock_src

# Copy the requirements file and install Python dependencies
COPY ./requirements.txt /app/
RUN pip install --default-timeout=100 -r requirements.txt

# Install additional dependencies for DiffDock
RUN pip install PyYAML scipy "networkx[default]" biopython rdkit-pypi e3nn spyrmsd pandas biopandas
RUN pip install "fair-esm[esmfold]"
RUN pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
RUN pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'

# Copy the application files into the container
COPY ./diffdock /app/diffdock

ENTRYPOINT ["uvicorn", "diffdock.api:app"]
