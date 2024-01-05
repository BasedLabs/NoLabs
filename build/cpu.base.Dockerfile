# Use an official PyTorch runtime as a parent image
FROM cnstark/pytorch:2.0.1-py3.9.17-ubuntu20.04

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update && \
    apt-get install -y cmake wget software-properties-common python3.8 python3.8-distutils curl libxrender1 libxext6 && \
    wget https://bootstrap.pypa.io/get-pip.py && \
    python3.8 get-pip.py

RUN pip install --upgrade cmake

RUN apt -y install default-jdk
RUN apt -y install default-jre

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
RUN apt-get install wget

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh
RUN conda --version
RUN conda create --name nolabs python=3.8