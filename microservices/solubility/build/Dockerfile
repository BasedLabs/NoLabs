# Use an official PyTorch runtime as a parent image
FROM cnstark/pytorch:2.0.1-py3.10.11-cuda11.8.0-ubuntu22.04

WORKDIR /app/solubility
WORKDIR /app
COPY requirements.txt /app
RUN pip install -r requirements.txt
COPY solubility/__init__.py /app/solubility
COPY solubility/predownload_models.py /app/solubility
COPY solubility/settings.ini /app/solubility
COPY solubility/settings.py /app/solubility
RUN python -m solubility.predownload_models
COPY . /app

ENTRYPOINT ["uvicorn", "solubility.api:app", "--workers", "2"]
