# Use a base image with CUDA support that matches the required PyTorch and CUDA versions
FROM pytorch/pytorch:1.11.0-cuda11.3-cudnn8-runtime

# Install necessary tools and libraries
RUN apt-get update && apt-get install -y wget git build-essential

# Set a working directory inside the container
WORKDIR /app

# Copy the requirements file and install Python dependencies
COPY ./requirements.txt /app/
# Ensure pip is up-to-date
RUN pip install --upgrade pip
# Install dependencies from requirements.txt
RUN pip install --default-timeout=100 -r requirements.txt

# Clone the DiffDock repository and checkout the specified commit
RUN git clone https://github.com/gcorso/DiffDock.git && \
    cd DiffDock && \
    git checkout a6c5275

# Install Python dependencies
RUN pip install --upgrade pip && \
    pip install pyg==0.7.1 \
                        pyyaml==6.0 \
                        scipy==1.7.3 \
                        networkx==2.6.3 \
                        biopython==1.79 \
                        rdkit-pypi==2022.03.5 \
                        e3nn==0.5.0 \
                        spyrmsd==0.5.2 \
                        pandas==1.5.3 \
                        biopandas==0.4.1

# Install PyTorch Geometric and its dependencies
RUN pip install git+https://github.com/pyg-team/pytorch_geometric.git
RUN pip install --no-index torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-1.11.0+cu113.html

# Clone the ESM repository, check out the specified commit, and install
RUN if [ ! -d "DiffDock/esm" ]; then \
        cd DiffDock && \
        git clone https://github.com/facebookresearch/esm && \
        cd esm && \
        git checkout ca8a710 && \
        pip install -e .; \
    fi

# Set up environment variables for ESM model weights
ENV HOME=/app/DiffDock/esm/model_weights

# Create directories for model weights and cache
RUN mkdir -p $HOME/.cache/torch/hub/checkpoints

# Download ESM model weights
RUN wget -P $HOME/.cache/torch/hub/checkpoints/ https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt && \
    wget -P $HOME/.cache/torch/hub/checkpoints/ https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50D-contact-regression.pt

RUN wget https://sourceforge.net/projects/smina/files/smina.static/download -O /app/DiffDock/smina && \
    chmod +x /app/DiffDock/smina

# Download and install gnina
RUN wget https://github.com/gnina/gnina/releases/download/v1.0.3/gnina -O /app/DiffDock/gnina && \
    chmod +x /app/DiffDock/gnina

# Copy the application files into the container
COPY ./diffdock /app/diffdock

ENTRYPOINT ["uvicorn", "diffdock.api:app"]