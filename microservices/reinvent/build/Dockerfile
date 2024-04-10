#FROM nvcr.io/nvidia/cuda:11.3.0-cudnn8-runtime-ubuntu20.04
FROM nvidia/cuda:12.3.2-runtime-ubuntu22.04

RUN apt-get update

RUN apt-get install -y wget libgomp1 unzip libxrender1 libxext6 git && rm -rf /var/lib/apt/lists/*

RUN wget -q \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /var/conda\
    && rm -f Miniconda3-latest-Linux-x86_64.sh

ENV PATH /var/conda/bin:$PATH

RUN conda --version

RUN conda install conda-libmamba-solver
RUN conda config --set solver libmamba
WORKDIR /app
RUN wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
RUN mv vina_1.2.5_linux_x86_64 /app/vina
RUN chmod +x /app/vina

# Install reinvent
WORKDIR /app/REINVENT4
RUN git clone https://github.com/MolecularAI/REINVENT4.git . && git checkout b88b2ba1a3e7d81480dc80ad403834084982628b
COPY install_reinvent.sh /app/REINVENT4
RUN sed -i 's/\r$//' install_reinvent.sh && chmod +x install_reinvent.sh
RUN apt-get update
RUN apt-get install build-essential -y
RUN ./install_reinvent.sh

# Install dockstream
WORKDIR /app/DockStream
RUN git clone https://github.com/BasedLabs/DockStream.git .
RUN conda env create -f /app/DockStream/environment.yml

WORKDIR /app

COPY reinvent /app/reinvent
RUN pip install -r reinvent/requirements.txt

COPY reinvent_configs /app/reinvent/reinvent_configs

RUN sed -i 's/\r$//' /app/reinvent/reinvent_configs/prepare_pdbqt.sh && chmod +x /app/reinvent/reinvent_configs/prepare_pdbqt.sh
RUN sed -i 's/\r$//' /app/reinvent/reinvent_configs/start_reinforcement_learning.sh && chmod +x /app/reinvent/reinvent_configs/start_reinforcement_learning.sh
RUN sed -i 's/\r$//' /app/reinvent/reinvent_configs/start_sampling.sh && chmod +x /app/reinvent/reinvent_configs/start_sampling.sh

ENTRYPOINT ["uvicorn", "reinvent.api:app"]