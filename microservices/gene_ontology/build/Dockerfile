# Use an official PyTorch runtime as a parent image
FROM cnstark/pytorch:2.0.1-py3.10.11-cuda11.8.0-ubuntu22.04

WORKDIR /app/gene_ontology
WORKDIR /app
COPY requirements.txt /app
RUN pip install -r requirements.txt
COPY gene_ontology/__init__.py /app/gene_ontology
COPY gene_ontology/predownload_models.py /app/gene_ontology
COPY gene_ontology/settings.ini /app/gene_ontology
COPY gene_ontology/settings.py /app/gene_ontology
RUN python -m gene_ontology.predownload_models
COPY . /app

ENTRYPOINT ["uvicorn", "gene_ontology.api:app"]
