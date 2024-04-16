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

NoLabs is an open source biolab that lets you run experiments with the latest state-of-the-art models for bio research.

The goal of the project is to accelerate bio research by making inference models easy to use for everyone. We are
currently supporting protein biolab (predicting useful protein properties such as solubility, localisation, gene
ontology, folding, etc.), drug discovery biolab (construct ligands and test binding to target proteins) and small molecules design biolab (design small molecules given a protein target and check drug-likeness and binding affinity).

We are working on expanding both and adding a cell biolab and genetic biolab, and we will appreciate your support and
contributions.

Let's accelerate bio research!

<img src="media/NoLabs_Architecture.png" width="100%">

# Features

**Bio Buddy - drug discovery co-pilot:**

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

To enable biobuddy run this command when starting nolabs:

```shell
$ ENABLE_BIOBUDDY=true docker compose up nolabs
```

And also start the biobuddy microservice:

```shell
$ OPENAI_API_KEY=your_openai_api_key TAVILY_API_KEY=your_tavily_api_key docker compose up biobuddy
```
Nolabs is running on GPT4 for the best performance. You can adjust the model you use in `microservices/biobuddy/biobuddy/services.py`

You can ignore OPENAI_API_KEY warnings when running other services using docker compose.

**Drug discovery lab:**

- Drug-target interaction prediction, high throughput virtual screening (HTVS) based on:
    - [DiffDock](https://github.com/gcorso/DiffDock)
    - [uMol](https://github.com/patrickbryant1/Umol)
- Automatic pocket prediction via [P2Rank](https://github.com/rdk/p2rank)
- Automatic MSA generation via [HH-suite3](https://github.com/soedinglab/hh-suite)

<br>
<img src="media/Docking.gif" width="100%">

**Protein lab:**

- Prediction of subcellular localisation via
  fine-tuned [ritakurban/ESM_protein_localization](https://huggingface.co/ritakurban/ESM_protein_localization) model (to
  be updated with a better model)
- Prediction of folded structure via [facebook/esmfold_v1](https://huggingface.co/facebook/esmfold_v1)
- Gene ontology prediction for 200 most popular gene ontologies
- Protein solubility prediction

<br>
<img src="media/localisation.gif" width="100%">

**Protein design Lab:**

- Protein generation via [RFDiffusion](https://github.com/RosettaCommons/RFdiffusion)

<br>
<img src="media/protein_design.gif" width="100%">

**Conformations Lab:**

- Conformations via [OpenMM](https://github.com/openmm/openmm) and [GROMACS](https://github.com/gromacs/gromacs)

**Small molecules design lab:**

- Small molecules design using a protein target with drug-likeness scoring component [REINVENT4](https://github.com/MolecularAI/REINVENT4)

Specify the search space (location) where designed molecule would bind relative to protein target. Then run reinforcement learning to generate new molecules in specified binding region.

WARNING: Reinforcement learning process might take a long time (with 128 molecules per 1 epoch and 50 epochs it could take a day)

<br>
<img src="media/small_molecules_design.png" width="100%">

# Starting

```bash
# Clone this project
$ git clone https://github.com/BasedLabs/nolabs
$ cd nolabs
```

Generate a new token for docker registry
https://github.com/settings/tokens/new
Select 'read:packages'

```bash
$ docker login ghcr.io -u username -p ghp_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

```bash
$ docker compose up
```

OR if you want to run a single feature

```bash
$ docker compose -up nolabs [gene_ontology|localisation|protein_design|solubility|conformations]
```

Server will be available on http://localhost:9000

# APIs

We provide individual Docker containers backed by FastAPI for each feature, which are available in the `/microservices`
folder. You can use them individually as APIs.

For example, to run the `esmfold` service, you can use Docker Compose:

```bash
$ docker compose up esmfold
```

Once the service is up, you can make a POST request to perform a task, such as predicting a protein's folded structure.
Here's a simple Python example:

```python
import requests

# Define the API endpoint
url = 'http://127.0.0.1:5736/run-folding'

# Specify the protein sequence in the request body
data = {
    'protein_sequence': 'YOUR_PROTEIN_SEQUENCE_HERE'
}

# Make the POST request and get the response
response = requests.post(url, json=data)

# Extract the PDB content from the response
pdb_content = response.json().get('pdb_content', '')

print(pdb_content)
```

This Python script makes a POST request to the esmfold microservice with a protein sequence and prints the predicted PDB
content.

## Running services on a separate machine

Since we provide individual Docker containers backed by FastAPI for each feature, available in the `/microservices`
folder, you can run them on separate machines. This setup is particularly useful if you're developing on a computer
without GPU support but have access to a VM with a GPU for tasks like folding, docking, etc.

For instance, to run the `diffdock` service, use Docker Compose on the VM or computer equipped with a GPU.

On your server/VM/computer with a GPU, run:

```bash
$ docker compose up diffdock
```

Once the service is up, you can check that you can access it from your computer by navigating to http://<
gpu_machine_ip>:5737/docs

If everything is correct, you should see the FastAPI page with diffdock's API surface like this:

<img src="media/Diffdock_fastapi.png">

Next, update the nolabs/infrastructure/settings.ini file on your primary machine to include the IP address of the
service (replace 127.0.0.1 with your GPU machine's IP):

```ini
...
p2rank = http://127.0.0.1:5731
esmfold = http://127.0.0.1:5736
esmfold_light = http://127.0.0.1:5733
msa_light = http://127.0.0.1:5734
umol = http://127.0.0.1:5735
diffdock = http://127.0.0.1:5737 -> http://74.82.28.227:5737
...
```

And now you are ready to use this service hosted on a separate machine!

## Supported microservices list

### 1) Protein design docker API

Model: [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion)

RFdiffusion is an open source method for structure generation, with or without conditional information (a motif, target
etc).

```shell
docker compose up protein_design
```

Swagger UI will be available on http://localhost:5789/docs

or
[install as a python package](microservices/protein_design/client/README.md)

### 2) ESMFold docker API

Model: [ESMFold](https://github.com/facebookresearch/esm) - Evolutionary Scale Modeling

```shell
docker compose up esmfold
```

Swagger UI will be available on http://localhost:5736/docs

or
[install as a python package](microservices/esmfold/client/README.md)

### 3) ESMAtlas docker API

Model: [ESMAtlas](https://esmatlas.com/about)

```shell
docker compose up esmfold_light
```

Swagger UI will be available on http://localhost:5733/docs

or
[install as a python package](microservices/esmfold_light/client/README.md)

### 4) Protein function prediction docker API

Model: [Hugging Face](https://huggingface.co/thomasshelby/go_prediction/resolve/main/go_model_150M.pth)

```shell
docker compose up gene_ontology
```

Swagger UI will be available on http://localhost:5788/docs

or
[install as a python package](microservices/gene_ontology/client/README.md)

### 5) Protein localisation prediction docker API

Model: [Hugging Face](https://huggingface.co/ritakurban/ESM_protein_localization)

```shell
docker compose up localisation
```

Swagger UI will be available on http://localhost:5787/docs

or
[install as a python package](microservices/localisation/client/README.md)

### 6) Protein binding site prediction docker API

Model: [p2rank](https://github.com/rdk/p2rank)

```shell
docker compose up p2rank
```

Swagger UI will be available on http://localhost:5731/docs

or
[install as a python package](microservices/p2rank/client/README.md)

### 7) Protein solubility prediction docker API

Model: [Hugging Face](https://huggingface.co/thomasshelby/solubility_model/resolve/main/solubility_model.pth)

```shell
docker compose up solubility
```

Swagger UI will be available on http://localhost:5786/docs

or
[Install as python package](microservices/solubility/client/README.md)

### 8) Protein-ligand structure prediction docker API

Model: [UMol](https://github.com/patrickbryant1/Umol)

```shell
docker compose up umol
```

Swagger UI will be available on http://localhost:5735/docs

or
[Install as python package](microservices/umol/client/README.md)

### 9) RoseTTAFold docker API

Model: [RoseTTAFold](https://github.com/RosettaCommons/RoseTTAFold)

```shell
docker compose up rosettafold
```

Swagger UI will be available on http://localhost:5738/docs

or
[Install as python package](microservices/rosettafold/client/README.md)

WARNING: To use Rosettafold you must change the <b>volumes '.'</b>
to point to the specified folders.

### 10) REINVENT4 Reinforcement Learning on a Protein receptor API

Model: [REINVENT4](https://github.com/MolecularAI/REINVENT4)

Misc: DockStream, QED, AutoDock Vina

```shell
docker compose up reinvent
```

Swagger UI will be available on http://localhost:5790/docs

or
[Install as python package](microservices/reinvent/client/README.md)

WARNING: Do not change the number of guvicorn workers (1), this will lead to microservice issues.

# Technologies

The following tools were used in this project:

- [Pytorch](https://pytorch.org/)
- [Jax](https://jax.readthedocs.io/en/latest/index.html)
- [Transformers](https://huggingface.co/transformers)
- [FastAPI](https://pypi.org/project/Flask/)
- [Docker](https://www.docker.com/)
- [Vue.js](https://vuejs.org/)
- [AutoDock Vina](https://vina.scripps.edu/)

## Requirements ##

**[Recommended for laptops]** If you are using a laptop, use ```--test``` argument (no need to have a lot of compute):

- RAM > 16GB
- [Optional] GPU memory >= 16GB (REALLY speeds up the inference)

**[Recommended for powerful workstations]** Else, if you want to host everything on your machine and have faster
inference (also a requirement for folding sequences > 400 amino acids in length):

- RAM > 30GB
- [Optional] GPU memory >= 40GB (REALLY speeds up the inference)

Made by <a href="https://github.com/jaktenstid" target="_blank">Igor</a>
and <a href="https://github.com/timurishmuratov7" target="_blank">Tim</a>

&#xa0;

<a href="#top">Back to top</a>
