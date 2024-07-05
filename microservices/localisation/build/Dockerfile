# Use an official PyTorch runtime as a parent image
FROM cnstark/pytorch:2.0.1-py3.10.11-cuda11.8.0-ubuntu22.04

WORKDIR /app/solubility
WORKDIR /app
COPY requirements.txt /app
RUN pip install -r requirements.txt
COPY . /app
RUN python -m localisation.predownload_model

ENTRYPOINT ["uvicorn", "localisation.api:app", "--workers", "2"]
