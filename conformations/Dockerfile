FROM continuumio/miniconda3:latest

WORKDIR /app
COPY gromacs.tar.gz /app

RUN apt-get -y update && \
    apt-get install -y cmake wget build-essential

RUN pip install --upgrade cmake

RUN ls -l
RUN tar xfz ./gromacs.tar.gz
WORKDIR /app/gromacs-2023.3
RUN mkdir build
WORKDIR /app/gromacs-2023.3/build
RUN cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
RUN make
RUN make check
RUN make install

ADD . /app

WORKDIR /app
RUN conda init zsh
RUN conda create --name conformations
RUN conda install -n conformations -c conda-forge -y --file requirements.txt

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "conformations", "python", "-u", "conformations_microservice.py", "--host=0.0.0.0", "--port=5785"]