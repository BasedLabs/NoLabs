<div align="center" id="top"> 
  <img src="NoLabs logo.png" alt="NoLabs" />

  &#xa0;

  <!-- <a href="https://nolabs.netlify.app">Demo</a> -->
</div>

<h1 align="center">NoLabs</h1>
<h2 align="center">Open source bioengine</h2>

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

NoLabs is an open source bioengine with support of web visualisation and hosting.

The goal of the project is to accelerate the Research stage of biotechnological processes via leveraging solutions for prediction and simulation.

## Features ##

1) Prediction of subcellular localisation via fine-tuned [ritakurban/ESM_protein_localization](https://huggingface.co/ritakurban/ESM_protein_localization) model (to be updated with a better model)

2) Prediction of folded structure via [facebook/esmfold_v1](https://huggingface.co/facebook/esmfold_v1) \

3) Inference of multiple proteins with saving to csv file

4) Docker containerisation for easy hosting

## Starting ##

```bash
# Clone this project
$ git clone https://github.com/BasedLabs/nolabs

# Access
$ cd nolabs

# Build the docker image
$ docker build -t nolabs .

# Run the image and expose the 5002 port
$ docker run -p 5002:5000 nolabs

# The website will be available on <http://localhost:5002>
```

## Technologies ##

The following tools were used in this project:

- [Pytorch](https://pytorch.org/)
- [Transformers](https://huggingface.co/transformers)
- [Flask](https://pypi.org/project/Flask/)
- [Docker](https://www.docker.com/)

## Requirements ##

- RAM > 32GB
- GPU memory >= 20GB (optional, but REALLY speeds up the folding inference)


## :memo: License ##

This project is under license from MIT. For more details, see the [LICENSE](LICENSE.md) file.


Made by <a href="https://github.com/jaktenstid" target="_blank">Igor</a> and <a href="https://github.com/timurishmuratov7" target="_blank">Tim</a>

&#xa0;

<a href="#top">Back to top</a>
