<meta>
<div align="center" id="top">
  <img src="media/NoLabs logo.png" alt="NoLabs" />
</div>

<h1 align="center">NoLabs</h1>
<h2 align="center">Open source biolab</h2>

<p align="center">
  <img alt="Github top language" src="https://img.shields.io/github/languages/top/BasedLabs/nolabs?color=56BEB8">
  <img alt="Github language count" src="https://img.shields.io/github/languages/count/BasedLabs/nolabs?color=56BEB8">
  <img alt="Repository size" src="https://img.shields.io/github/repo-size/BasedLabs/nolabs?color=56BEB8">
  <img alt="License" src="https://img.shields.io/github/license/BasedLabs/nolabs?color=56BEB8">
</p>

# Contents

* [About](#about)
* [Features](#features)
* [Starting](#starting)
* [Microservices](#microservices)
* [Technologies](#technologies)

# About

NoLabs is an open source biolab that lets you run experiments with the latest state-of-the-art models, bioinformatics tools and scalable no-code workflow engine for bio research.

<img src="media/NoLabs_Architecture.png" width="100%">

# Features

**Workflow Engine**

<br>
<img src="media/Workflow_UI_demo.png" width="100%">

- Create workflows combining different models and data
- Schedule jobs and observe results for big data processing
- Adjust input parameters for particular jobs

# Starting

⚠️ **Warning:** For macOS users, there are known issues with running Docker Compose properly for certain setups. Please use COLIMA: https://github.com/abiosoft/colima

```bash
# Clone this project
$ git clone https://github.com/BasedLabs/nolabs
$ cd nolabs
# Create .env files (you will be able to adjust them)
$ chmod +x scripts/gen-envs.sh
$ make gen-envs
```

Generate a new token for docker registry
https://github.com/settings/tokens/new
Select 'read:packages'

```bash
$ docker login ghcr.io -u username -p ghp_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

If you want to run a single feature **(recommended)**

```bash
# use gpu or no-gpu profile
$ docker compose up nolabs-frontend nolabs-api nolabs-worker mongo redis
# mongo, redis and worker are required
$ docker compose up esmfold-light-gpu
$ docker compose up diffdock-gpu
...
```

OR if you want to run everything on one machine:

```bash
$ docker compose --profile gpu up
# or
$ docker compose --profile cpu up
```
Server will be available on http://localhost:9000

These commands will start

- NoLabs frontend (appication UI server)
- NoLabs worker (fault-tolerant distributed workflow scheduler)
- NoLabs API backend (API backend for UI)

Since workflow is backend by Celery and Redis, <b>Redis instance must be accessible from all workers!</b>

# Scalability

Workflow engine is scalable thanks to Celery and Redis, so you can run multiple instances of, let's say, Diffdock and jobs will be automatically distributed on multiple machines. By the default each worker processes one job at a time.

# Settings

You can find `.env` files in `nolabs/infrastructure/.env` and under `microservices` folder. You can adjust settings as you wish.
To regenerate `.env` files from templates run 

```shell
make gen-envs
```

# APIs

We provide individual Docker containers backed by FastAPI for each feature, which are available in the `/microservices`
folder. You can use them individually as APIs.

For example, to run the `esmfold` service, you can use Docker Compose:

```bash
$ docker compose up esmfold-api
# or
$ docker compose up esmfold-gpu-api
```

Once the service is up, you can make a POST request to perform a task, such as predicting a protein's folded structure.
Here's a simple Python example:

```python
import requests

# Define the API endpoint
url = 'http://127.0.0.1:5736/inference'

# Specify the protein sequence in the request body
data = {
    'fasta_sequence': ''
}

# Make the POST request and get the response
response = requests.post(url, json=data)

# Extract the PDB content from the response
pdb_content = response.json().get('pdb_content', '')

print(pdb_content)
```

This Python script makes a POST request to the esmfold microservice with a protein sequence and prints the predicted PDB
content.

## Running workers/API on a separate machine

Since we provide individual Docker containers backed by Celery\Fastapi for each feature, available in the `microservices`
folder, you can run them on separate machines. This setup is particularly useful if you're developing on a computer
without GPU support but have access to a VM with a GPU for tasks like folding, docking, etc.

For instance, to run the `diffdock` service, use Docker Compose on the VM or computer equipped with a GPU.

On machine with a GPU do the following steps:

`   1) Adjust `microservices/diffdock/service/.env` file with `REDIS_URL` pointing to your Redis deployment
2) Run following command
```bash
$ docker compose up diffdock-gpu
```

You should see celery start log message.
And now you are ready to use this service hosted on a separate machine!

In case you want to start just an API (there is no integration in NoLabs, however you can use service separately):

```bash
$ docker compose up diffdock-api-gpu
```

Once the service is up, you can check that you can access it from your computer by navigating to `http://<
gpu_machine_ip>:5737/docs`

## Bio Buddy - drug discovery co-pilot

BioBuddy is a drug discovery copilot that supports:

- Downloading data from [ChemBL](https://www.ebi.ac.uk/chembl/)
- Downloading data from [RcsbPDB](https://www.rcsb.org/)
- Questions about drug discovery process, targets, chemical components etc
- Writing review reports based on published papers

For example, you can ask
- "Can you pull me some latest approved drugs?"
- "Can you download me 1000 rhodopsins?"
- "How does an aspirin molecule look like?" and it will do this and answer other questions.

<img src="media/Biobuddy_pic.png" width="100%">

To enable biobuddy:

1) Adjust `NOLABS_ENABLE_BIOBUDDY` environment variable in `nolabs/infrastructure/.env`
2) Adjust `OPENAI_API_KEY` and `TAVILY_API_KEY` in `microservices/biobuddy/biobuddy/.env`
3) Start docker compose
```shell
$ docker compose up --profile gpu biobuddy nolabs-frontend nolabs-worker nolabs-api mongo redis
```

Nolabs is running on GPT4 for the best performance. You can adjust the model you use in `microservices/biobuddy/biobuddy/services.py`


## How-to run workers

### 1) RFdiffusion (protein design)

Model: [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion)

RFdiffusion is an open source method for structure generation, with or without conditional information (a motif, target
etc).

```shell
make download-rfdiffusion-weights
docker compose up rfdiffusion-gpu
```

### 2) ESMFold (folding)

Model: [ESMFold](https://github.com/facebookresearch/esm) - Evolutionary Scale Modeling

```shell
make download-esmfold-weights
docker compose up esmfold-gpu
```

### 3) ESMAtlas (folding)

Model: [ESMAtlas](https://esmatlas.com/about)

```shell
docker compose up esmfold-light
```

### 4) Diffdock (protein-ligand binding prediction)

Model: [DiffDock](https://github.com/gcorso/DiffDock)

```shell
make download-diffdock-weights
docker compose up diffdock-gpu
```

### 4) Proteinmpnn (design fasta from pdb)

Model: [Proteinmpnn](https://github.com/dauparas/ProteinMPNN)

```shell
docker compose up proteinmpnn-gpu
```

## Requirements ##

**[Recommended for laptops]**

- RAM > 16GB
- [Optional] GPU memory >= 16GB (REALLY speeds up the inference)

**[Recommended for powerful workstations]** Else, if you want to host everything on your machine and have faster
inference (also a requirement for folding sequences > 400 amino acids in length):

- RAM > 30GB
- [Optional] GPU memory >= 40GB (REALLY speeds up the inference)

Made by <a href="https://github.com/jaktenstid" target="_blank">Igor Bruev</a>
and <a href="https://github.com/timurishmuratov7" target="_blank">Tim Ishmuratov</a>

&#xa0;

<a href="#top">Back to top</a>
