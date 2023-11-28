# Use an official PyTorch runtime as a parent image
FROM pytorch/pytorch:2.1.0-cuda12.1-cudnn8-devel

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && \
    apt-get install -y cmake wget software-properties-common python3.8 python3.8-distutils curl libxrender1 libxext6 && \
    wget https://bootstrap.pypa.io/get-pip.py && \
    python3.8 get-pip.py

# Continue with your setup...
WORKDIR /app
ADD . /app

RUN tar xfz gromacs.tar.gz
WORKDIR /app/gromacs-2023.3
RUN mkdir build
RUN cd build
RUN cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
RUN make
RUN make check
RUN sudo make install
RUN source /usr/local/gromacs/bin/GMXRC

WORKDIR /app

# Install Node.js 18
RUN apt-get install -y ca-certificates curl gnupg
RUN mkdir -p /etc/apt/keyrings
RUN curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg
ARG NODE_MAJOR=18
RUN echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_$NODE_MAJOR.x nodistro main" | tee /etc/apt/sources.list.d/nodesource.list
RUN apt-get update
RUN apt-get install nodejs -y

RUN apt-get install libstdc++6

RUN python3.8 -m pip install --upgrade pip && \
    python3.8 -m pip install torch==1.13.1 torchvision torchaudio && \
    python3.8 -m pip install torch-geometric && \
    python3.8 -m pip install gcc=12.1.0 openmm pdbfixer parmed && \
    python3.8 -m pip install --no-index pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-1.13.0+cu116.html && \
    python3.8 -m pip install -r requirements.txt && \

WORKDIR /app/src/server/frontend
RUN rm -rf node_modules
RUN npm install
WORKDIR /app
EXPOSE 5000
EXPOSE 5137
ENV PYTHONPATH=/app
RUN chmod +x ./dockerfile_entrypoint.sh
ENTRYPOINT ["sh", "./dockerfile_entrypoint.sh"]