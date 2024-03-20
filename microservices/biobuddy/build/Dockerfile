# Use a PyTorch base image with CUDA support
FROM python:3.8

# Install Python 3.8 and curl
RUN apt-get update && apt-get install -y curl python3-pip

# Set a working directory
WORKDIR /app

# Set the OPENAI_API_KEY environment variable
ENV OPENAI_API_KEY=abc123

# Copy only the requirements.txt first to leverage Docker cache
COPY ./requirements.txt /app/

# Install dependencies
RUN python3.8 -m pip install --upgrade pip
RUN python3.8 -m pip install --default-timeout=100 -r requirements.txt

# Now copy the rest of the application files into the container
COPY ./biobuddy /app/biobuddy

# Set the command to run the FastAPI app
ENTRYPOINT ["uvicorn", "biobuddy.api:app"]