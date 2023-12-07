# Use an official PyTorch runtime as a parent image
FROM ghcr.io/basedlabs/nolabs-cpu:latest

# Continue with your setup...
ADD . /app
WORKDIR /app

RUN python3.8 -m pip install --upgrade pip && \
    python3.8 -m pip install torch==1.13.1 torchvision torchaudio && \
    python3.8 -m pip install torch-geometric && \
    python3.8 -m pip install --no-index pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-1.13.0+cu116.html

RUN python3.8 -m pip install -r requirements.txt

WORKDIR /app/src/server/frontend
RUN rm -rf node_modules
RUN npm install
WORKDIR /app
EXPOSE 5000
EXPOSE 5137
ENV PYTHONPATH=/app
RUN chmod +x ./dockerfile_entrypoint.sh
ENTRYPOINT ["sh", "./dockerfile_entrypoint.sh"]