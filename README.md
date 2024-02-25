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

## About ##

NoLabs is an open source biolab which lets you run experiments with latest state of the art models for bio research.

The goal of the project is to accelerate bio research by making inference models easy to use for everyone. We are currenly supporting protein biolab (predicting useful protein properties such as solubility, localisation, gene ontology, folding etc.) and drug discovery biolab (construct ligands and test binding to target proteins). 

We are working on expanding both and adding a cell biolab and genetic biolab, and we will appreciate your support and contributions. 

Let's accelerate bio research!

<img src="media/NoLabs_Architecture.png" width="100%">


## Features ##

**Drug discovery lab (State of the art):**
- Drug-target interaction prediction, high throughput virtual screening (HTVS) based on:
  - [DiffDock](https://github.com/gcorso/DiffDock)
  - [uMol](https://github.com/patrickbryant1/Umol)
- Automatic pocket prediction via [P2Rank](https://github.com/rdk/p2rank)
- Automatic MSA generation via [HH-suite3](https://github.com/soedinglab/hh-suite)

<br>
<img src="media/Docking.gif" width="100%">

**Protein lab:**

- Prediction of subcellular localisation via fine-tuned [ritakurban/ESM_protein_localization](https://huggingface.co/ritakurban/ESM_protein_localization) model (to be updated with a better model)
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

## Starting ##

```bash
# Clone this project
$ git clone https://github.com/BasedLabs/nolabs
$ cd nolabs
```

```bash
$ docker compose up
```
OR if you want to run a single feature

```bash
$ docker compose -up nolabs [gene_ontology|localisation|protein_design|solubility|conformations]
```

Server will be available on http://localhost:9000

## APIs ##

We provide individual Docker containers backed by FastAPI for each feature, which are available in the `/microservices` folder. You can use them individually as APIs.

For example, to run the `esmfold` service, you can use Docker Compose:

```bash
$ docker compose up esmfold
```

Once the service is up, you can make a POST request to perform a task, such as predicting a protein's folded structure. Here's a simple Python example:

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
This Python script makes a POST request to the esmfold microservice with a protein sequence and prints the predicted PDB content.

## Running services on a separate machine 

Since we provide individual Docker containers backed by FastAPI for each feature, available in the `/microservices` folder, you can run them on separate machines. This setup is particularly useful if you're developing on a computer without GPU support but have access to a VM with a GPU for tasks like folding, docking, etc.

For instance, to run the `diffdock` service, use Docker Compose on the VM or computer equipped with a GPU. 

On your server/VM/computer with a GPU, run:

```bash
$ docker compose up diffdock
```

Once the service is up, you can check that you can access it from your computer by navigating to http://<gpu_machine_ip>:5737/docs

If everything is correct, you should see the FastAPI page with diffdock's API surface like this:

<img src="media/Diffdock_fastapi.png">

Next, update the nolabs/infrastructure/settings.ini file on your primary machine to include the IP address of the service (replace 127.0.0.1 with your GPU machine's IP):

```ini
...
p2rank=http://127.0.0.1:5731
esmfold=http://127.0.0.1:5736
esmfold_light=http://127.0.0.1:5733
msa_light=http://127.0.0.1:5734
umol=http://127.0.0.1:5735
diffdock=http://127.0.0.1:5737 -> http://74.82.28.227:5737
...
```

And now you are ready to use this service hosted on a separate machine!

## Technologies ##

The following tools were used in this project:

- [Pytorch](https://pytorch.org/)
- [Jax](https://jax.readthedocs.io/en/latest/index.html)
- [Transformers](https://huggingface.co/transformers)
- [FastAPI](https://pypi.org/project/Flask/)
- [Docker](https://www.docker.com/)
- [Vue.js](https://vuejs.org/)

## Requirements ##

**[Recommended for laptops]** If you are using a laptop, use ```--test``` argument (no need to have a lot of compute):
- RAM > 16GB
- [Optional] GPU memory >= 16GB (REALLY speeds up the inference)

**[Recommended for powerful workstations]** Else, if you want to host everything on your machine and have faster inference (also a requirement for folding sequences > 400 amino acids in length):
- RAM > 30GB
- [Optional] GPU memory >= 40GB (REALLY speeds up the inference)

Made by <a href="https://github.com/jaktenstid" target="_blank">Igor</a> and <a href="https://github.com/timurishmuratov7" target="_blank">Tim</a>

&#xa0;

<a href="#top">Back to top</a>
