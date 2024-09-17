# Use a lightweight Python base image with necessary tools
FROM python:3.9-slim

# Install necessary tools and libraries
RUN apt-get update && apt-get install -y \
    wget \
    git \
    build-essential \
    curl \
    unzip

# Set a working directory inside the container
WORKDIR /app

# Install Python dependencies from a requirements file (if any)
COPY ./requirements.txt /app/
RUN pip install --upgrade pip
RUN pip install --default-timeout=100 -r requirements.txt

# Install MMseqs2
RUN wget https://mmseqs.com/latest/mmseqs-linux.tar.gz && \
    tar -xzf mmseqs-linux.tar.gz && \
    mv mmseqs /usr/local/bin && \
    rm mmseqs-linux.tar.gz

# Set up MMseqs2 database directory
RUN mkdir -p /app/mmseqs2_databases

# Download MMseqs2 databases (you can customize which database you need)
RUN wget https://mmseqs.com/uniref/uniref90/uniref90_03_2021.tar.gz -P /app/mmseqs2_databases/ && \
    tar -xzf /app/mmseqs2_databases/uniref90_03_2021.tar.gz -C /app/mmseqs2_databases/ && \
    rm /app/mmseqs2_databases/uniref90_03_2021.tar.gz

# Make sure MMseqs2 is accessible in the environment
ENV PATH="/usr/local/bin/mmseqs:${PATH}"

# Add a script for handling MSA queries (can be a REST API or CLI tool)
COPY ./msa_service /app/msa

# Set entry point for running your service (e.g., uvicorn if using FastAPI, or python app.py for a Flask app)
ENTRYPOINT ["uvicorn", "msa.api:app"]
