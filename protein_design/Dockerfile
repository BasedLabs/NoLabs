FROM ghcr.io/basedlabs/rfdiffusion:latest

ENV HYDRA_FULL_ERROR=1
RUN apt-get -y install git

WORKDIR /app
RUN apt-get update
RUN apt-get -y install wget
ADD . /app
RUN mkdir RFdiffusion/models
RUN chmod +x RFdiffusion/scripts/download_models.sh
RUN RFdiffusion/scripts/download_models.sh /app/RFdiffusion/models
ADD requirements.txt /app
RUN pip install -r requirements.txt
RUN ls -l

ENTRYPOINT ["python3.9", "-u", "protein_design_microservice.py", "--host=0.0.0.0", "--port=5786"]