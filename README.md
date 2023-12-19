<div align="center" id="top"> 
  <img src="media/NoLabs logo.png" alt="NoLabs" />

  &#xa0;

  <!-- <a href="https://nolabs.netlify.app">Demo</a> -->
</div>

<h1 align="center">NoLabs</h1>
<h2 align="center">Open source biolab</h2>

<p align="center">
  <img alt="Github top language" src="https://img.shields.io/github/languages/top/BasedLabs/nolabs?color=56BEB8">

  <img alt="Github language count" src="https://img.shields.io/github/languages/count/BasedLabs/nolabs?color=56BEB8">

  <img alt="Repository size" src="https://img.shields.io/github/repo-size/BasedLabs/nolabs?color=56BEB8">

  <img alt="License" src="https://img.shields.io/github/license/BasedLabs/nolabs?color=56BEB8">

  <!-- <img alt="Github issues" src="https://img.shields.io/github/issues/BasedLabs/nolabs?color=56BEB8" /> -->

  <!-- <img alt="Github forks" src="https://img.shields.io/github/forks/BasedLabs/nolabs?color=56BEB8" /> -->

  <!-- <img alt="Github stars" src="https://img.shields.io/github/stars/BasedLabs/nolabs?color=56BEB8" /> -->
</p>

<!-- Status -->

<!-- <h4 align="center"> 
	🚧  NoLabs 🚀 Under construction...  🚧
</h4> 

<hr> -->

<p align="center">
  <a href="#dart-about">About</a> &#xa0; | &#xa0; 
  <a href="#sparkles-features">Features</a> &#xa0; | &#xa0;
  <a href="#rocket-technologies">Technologies</a> &#xa0; | &#xa0;
  <a href="#white_check_mark-requirements">Requirements</a> &#xa0; | &#xa0;
  <a href="#checkered_flag-starting">Starting</a> &#xa0; | &#xa0;
  <a href="#memo-license">License</a> &#xa0; | &#xa0;
  <a href="https://github.com/BasedLabs" target="_blank">Author</a>
</p>

<br>

## About ##

NoLabs is an open source biolab with support of web visualisation and hosting.

The goal of the project is to accelerate bio research via making inference models easy to use for everyone. We are currenly supporting protein biolab (predicting useful protein properties such as solubility, localisation, Gene ontology, folding etc.) and drug discovery biolab (construct ligands and test binding to target proteins). 

We are working on expanding both and adding a cell biolab and genetic biolab, and we will appreciate your support and contributions. Let's accelerate bio research!

<hr style="border:2px solid gray">
<img src="media/website-screenshot.jpg" width="100%">
<hr style="border:2px solid gray">

[<img src="media/amino-acid-lab.gif" width="100%">](media/NoLabs.mp4)
[<img src="media/drug-target.gif" width="100%">](media/NoLabs.mp4)

<hr style="border:2px solid gray">

## Features ##

**Protein biolab:**

1) Prediction of subcellular localisation via fine-tuned [ritakurban/ESM_protein_localization](https://huggingface.co/ritakurban/ESM_protein_localization) model (to be updated with a better model)

2) Prediction of folded structure via [facebook/esmfold_v1](https://huggingface.co/facebook/esmfold_v1)

3) Gene ontology prediction for 200 most popular gene ontologies

4) Protein solubility prediction

**Drug discovery biolab:**
1) Drug-target interaction prediction, hight throughput virtual screening (HTVS) based on [uMol](https://github.com/patrickbryant1/Umol)

Hosting:
1) Docker containerisation for easy hosting

## Starting ##

```bash
# Clone this project
$ git clone https://github.com/BasedLabs/nolabs

# Access
$ cd nolabs

# Build the docker image
$ docker build -t nolabs -f conformations.Dockerfile .
$ docker buil 

# Run the image and expose the 5000 port
$ docker run -p 5173:5173 -p 5000:5000 nolabs --test
# Run without --test if you want to run models on GPU
# check 'Requirements' section for more information

# The website will be available on <http://localhost:5173>
```

## Technologies ##

The following tools were used in this project:

- [Pytorch](https://pytorch.org/)
- [Jax](https://jax.readthedocs.io/en/latest/index.html)
- [Transformers](https://huggingface.co/transformers)
- [Flask](https://pypi.org/project/Flask/)
- [Docker](https://www.docker.com/)

## Requirements ##

**[Recommended for laptops]** If you are using a laptop, use ```--test``` argument (no need to have a lot of compute):
- RAM > 16GB
- [Optional] GPU memory >= 16GB (REALLY speeds up the inference)

**[Recommended for powerful workstations]** Else, if you want to host everything on your machine and have faster inference:
- RAM > 30GB
- [Optional] GPU memory >= 40GB (REALLY speeds up the inference)

## :memo: License ##

This project is under license from MIT. For more details, see the [LICENSE](LICENSE.md) file.


Made by <a href="https://github.com/jaktenstid" target="_blank">Igor</a> and <a href="https://github.com/timurishmuratov7" target="_blank">Tim</a>

&#xa0;

<a href="#top">Back to top</a>
