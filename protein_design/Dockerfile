FROM ghcr.io/basedlabs/rfdiffusion:latest

ENV HYDRA_FULL_ERROR=1

WORKDIR /app
ADD requirements.txt /app
RUN pip install -r requirements.txt
ADD . /app
RUN ls -l

ENTRYPOINT ["python3.9", "-u", "protein_design_microservice.py", "--host=0.0.0.0", "--port=5786"]