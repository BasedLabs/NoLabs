# Use a base image with CUDA support that matches the required PyTorch and CUDA versions
FROM pytorch/pytorch:1.11.0-cuda11.3-cudnn8-runtime

# Install necessary tools and libraries
RUN apt-get update && apt-get install -y wget git build-essential

# Set a working directory inside the container
WORKDIR /app

# Copy the requirements file and install Python dependencies
COPY ./requirements.txt /app/
# Ensure pip is up-to-date
RUN pip install --upgrade pip
# Install dependencies from requirements.txt
RUN pip install --default-timeout=100 -r requirements.txt

# Clone the ProtMPNN repository
RUN git clone https://github.com/dauparas/ProteinMPNN.git /app/ProtMPNN

# Copy the application files into the container
COPY ./prot_mpnn /app/prot_mpnn

# Expose the port (if needed)
EXPOSE 8000

# Set the entrypoint to start the FastAPI app
ENTRYPOINT ["uvicorn", "prot_mpnn.api:app", "--host", "0.0.0.0", "--port", "8000"]
