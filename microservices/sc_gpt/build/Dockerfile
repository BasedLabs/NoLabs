FROM python:3.10.12

# Install Python and curl
RUN apt-get update && apt-get install -y curl python3-pip gfortran build-essential libopenblas-dev
RUN apt-get install -y libhdf5-dev

# Set a working directory
WORKDIR /app

# Copy only the requirements.txt first to leverage Docker cache
COPY ./requirements.txt /app/

# Install dependencies
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install -r requirements.txt

# Download the model
COPY ./sc_gpt/download_models.py /app/sc_gpt/download_models.py
RUN python3 /app/sc_gpt/download_models.py

# Now copy the rest of the application files into the container
COPY ./sc_gpt /app/sc_gpt

# Set the command to run the FastAPI app
ENTRYPOINT ["uvicorn", "sc_gpt.api:app"]
