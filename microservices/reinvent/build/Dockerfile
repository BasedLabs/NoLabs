FROM nvcr.io/nvidia/cuda:11.3.0-cudnn8-runtime-ubuntu20.04
ARG ROSETTACOMMONS_CONDA_USERNAME
ARG ROSETTACOMMONS_CONDA_PASSWORD

RUN apt-get update

RUN apt-get install -y wget libgomp1 unzip && rm -rf /var/lib/apt/lists/*

RUN wget -q \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /var/conda\
    && rm -f Miniconda3-latest-Linux-x86_64.sh

ENV PATH /var/conda/bin:$PATH

RUN conda --version

RUN conda install conda-libmamba-solver
RUN conda config --set solver libmamba
WORKDIR /app
COPY reinvent/ReinventCommunity /app
RUN conda env create -f environment.yml
WORKDIR /app/temp
COPY install.sh /app/temp
RUN sed -i 's/\r$//' install.sh && chmod +x install.sh

RUN ./install.sh

RUN mv vina_1.2.5_linux_x86_64 /app
WORKDIR /app
COPY run.sh /app
WORKDIR /app/REINVENT4
WORKDIR /app/DockStream
WORKDIR /app
COPY /reinvent/REINVENT4 /app/REINVENT4
COPY /reinvent/DockStream /app/DockStream
RUN conda env create -f /app/DockStream/environment.yml
COPY install_reinvent.sh /app/REINVENT4
WORKDIR /app/REINVENT4
RUN sed -i 's/\r$//' install_reinvent.sh && chmod +x install_reinvent.sh
RUN ./install_reinvent.sh

RUN apt-get update
RUN apt-get install libxrender1 libxext6 -y
WORKDIR /app
RUN sed -i 's/\r$//' run.sh && chmod +x run.sh
COPY example.toml /app
COPY case_study_files /app/case_study_files
COPY 1smd.pdb /app
ENTRYPOINT ["./run.sh"]