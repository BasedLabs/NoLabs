# Use an official PyTorch runtime as a parent image
FROM pytorch/pytorch:2.1.0-cuda12.1-cudnn8-devel

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update && \
    apt-get install -y cmake wget software-properties-common python3.8 python3.8-distutils curl libxrender1 libxext6 && \
    wget https://bootstrap.pypa.io/get-pip.py && \
    python3.8 get-pip.py

RUN pip install --upgrade cmake

# Continue with your setup...
WORKDIR /app
ADD . /app

RUN tar xfz gromacs.tar.gz
WORKDIR /app/gromacs-2023.3
RUN mkdir build
WORKDIR /app/gromacs-2023.3/build
RUN cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
RUN make
RUN make check
RUN make install

WORKDIR /app

# Install Node.js 18
RUN apt-get install -y ca-certificates curl gnupg
RUN mkdir -p /etc/apt/keyrings
RUN curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg
ARG NODE_MAJOR=18
RUN echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_$NODE_MAJOR.x nodistro main" | tee /etc/apt/sources.list.d/nodesource.list
RUN apt-get update
RUN apt-get install nodejs -y

RUN apt-get install libstdc++6 -y