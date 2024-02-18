# Use an official PyTorch runtime as a parent image
FROM cnstark/pytorch:2.0.1-py3.10.11-cuda11.8.0-ubuntu22.04

RUN apt-get update -y
RUN apt-get install -y git

WORKDIR /microservice
COPY reinvent/ /microservice/reinvent
COPY requirements.txt /microservice
RUN ls -l
RUN pip install -r reinvent/REINVENT4/requirements-linux-64.lock
RUN pip install --no-deps reinvent/REINVENT4
RUN apt-get install -y libxrender1
RUN reinvent --help

RUN mkdir -p ~/miniconda3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
RUN bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
RUN rm -rf ~/miniconda3/miniconda.sh
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN conda env create -f reinvent/DockStream/environment.yml

RUN wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
RUN chmod +x vina_1.2.5_linux_x86_64
RUN ls -l
RUN ./vina_1.2.5_linux_x86_64

COPY app/ /microservice/app

RUN pip install -r requirements.txt

ENTRYPOINT ["uvicorn", "app.api:app", "--host=0.0.0.0"]
