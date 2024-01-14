FROM ghcr.io/basedlabs/conformations-base:latest

COPY requirements.txt /app
WORKDIR /app

RUN conda init zsh
RUN conda create --name conformations
RUN conda install -n conformations -c conda-forge -y --file requirements.txt
COPY . /app
WORKDIR /app

ENV PYTHONBUFFERED=1

# Append host and port to the end
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "conformations", "uvicorn", "conformations.api:app"]
