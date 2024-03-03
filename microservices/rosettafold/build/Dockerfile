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

RUN apt-get update
RUN apt-get -y install wget git
WORKDIR /RoseTTAFold
RUN git clone https://github.com/RosettaCommons/RoseTTAFold.git /RoseTTAFold
COPY ./run_e2e_ver.sh /RoseTTAFold

RUN conda install conda-libmamba-solver
RUN conda config --set solver libmamba
RUN conda env create -q -f RoseTTAFold-linux.yml
RUN conda env create -q -f folding-linux.yml

RUN conda config --add channels https://${ROSETTACOMMONS_CONDA_USERNAME}:${ROSETTACOMMONS_CONDA_PASSWORD}@conda.graylab.jhu.edu
#installing pyrosetta into a base image so it gets cached between builds
RUN conda install -n folding -y pyrosetta

RUN wget -q https://files.ipd.uw.edu/pub/RoseTTAFold/weights.tar.gz --show-progress
RUN tar xfz weights.tar.gz
##
RUN sed -i 's/\r$//' install_dependencies.sh && chmod +x install_dependencies.sh
RUN ./install_dependencies.sh

ENV PATH /RoseTTAFold:$PATH

COPY microservice /RoseTTAFold/microservice
RUN sed -i 's/\r$//' run_pyrosetta_ver.sh && chmod +x run_pyrosetta_ver.sh
RUN sed -i 's/\r$//' input_prep/make_msa.sh && chmod +x input_prep/make_msa.sh
RUN sed -i 's/\r$//' input_prep/make_ss.sh && chmod +x input_prep/make_ss.sh

COPY requirements.txt /RoseTTAFold
RUN pip install -r requirements.txt
ENTRYPOINT ["uvicorn", "microservice.api:app"]