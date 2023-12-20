FROM ghcr.io/basedlabs/rfdiffusion:latest

WORKDIR /app
ADD . /app
RUN pip install -r requirements.txt

ENTRYPOINT ["python", "-u", "protein_design_microservice.py", "--host=0.0.0.0", "--port=5786"]