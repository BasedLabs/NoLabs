# Use an official PyTorch runtime as a parent image
FROM pytorch/pytorch:2.1.0-cuda12.1-cudnn8-devel

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && \
    apt-get install -y wget software-properties-common python3.8 python3.8-distutils curl libxrender1 libxext6 && \
    wget https://bootstrap.pypa.io/get-pip.py && \
    python3.8 get-pip.py

# Install Node.js 18
RUN curl -sL https://deb.nodesource.com/setup_16.x | bash -
RUN apt-get install -y nodejs

# Continue with your setup...
WORKDIR /app
ADD . /app

RUN python3.8 -m pip install --upgrade pip && \
    python3.8 -m pip install torch==1.13.1 torchvision torchaudio && \
    python3.8 -m pip install torch-geometric && \
    python3.8 -m pip install --no-index pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-1.13.0+cu116.html && \
    python3.8 -m pip install -r requirements.txt

WORKDIR /app/src/server/frontend
RUN rm -rf node_modules
RUN npm install
WORKDIR /app
EXPOSE 5000
EXPOSE 5137
ENV PYTHONPATH=/app
RUN chmod +x ./dockerfile_entrypoint.sh
ENTRYPOINT ["sh", "./dockerfile_entrypoint.sh"]