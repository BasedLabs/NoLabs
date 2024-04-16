# Use an official PyTorch runtime as a parent image
FROM malevichai/app:python-torch_v0.1

RUN mkdir ./apps
COPY requirements.txt ./apps
RUN pip install -r ./apps/requirements.txt
COPY . ./apps
RUN (cd ./apps && python -m localisation.predownload_model)

ENV PYTHONPATH "${PYTHONPATH}:/julius/apps"