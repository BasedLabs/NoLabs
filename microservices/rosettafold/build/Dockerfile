FROM cnstark/pytorch:2.0.1-py3.9.17-cuda11.8.0-devel-ubuntu20.04

RUN apt-get update

RUN apt-get install -y wget libgomp1 unzip tar && rm -rf /var/lib/apt/lists/*

RUN wget -q \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /var/conda\
    && rm -f Miniconda3-latest-Linux-x86_64.sh

ENV PATH /var/conda/bin:$PATH

RUN conda --version

RUN apt-get update
RUN apt-get -y install wget git vim
WORKDIR /RoseTTAFold
RUN git clone https://github.com/uw-ipd/RoseTTAFold2.git /RoseTTAFold && git checkout 703d56e
COPY RF2-linux.yml /RoseTTAFold
RUN conda env create -f RF2-linux.yml
COPY install_transformers.sh /RoseTTAFold
RUN sed -i 's/\r$//' install_transformers.sh && chmod +x install_transformers.sh
RUN ./install_transformers.sh

WORKDIR /RoseTTAFold/network
RUN wget https://files.ipd.uw.edu/dimaio/RF2_apr23.tgz
RUN tar xvfz RF2_apr23.tgz
WORKDIR /RoseTTAFold

COPY run_RF2.sh /RoseTTAFold
RUN sed -i 's/\r$//' run_RF2.sh && chmod +x run_RF2.sh
ENV PATH /RoseTTAFold:$PATH

COPY run_RF2_msa.sh /RoseTTAFold
RUN sed -i 's/\r$//' run_RF2_msa.sh && chmod +x run_RF2_msa.sh
ENV PATH /RoseTTAFold:$PATH

COPY microservice /RoseTTAFold/microservice
COPY requirements.txt /RoseTTAFold
RUN pip install -r requirements.txt

ENTRYPOINT ["uvicorn", "microservice.api:app"]
