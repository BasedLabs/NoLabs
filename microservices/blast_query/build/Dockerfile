FROM python:3.11

# Install Python 3.11 and curl
RUN apt-get update && apt-get install -y curl python3-pip

# Set a working directory
WORKDIR /app

# Copy only the requirements.txt first to leverage Docker cache
COPY ./requirements.txt /app/

# Install dependencies
RUN pip install --upgrade pip
RUN pip install --default-timeout=100 -r requirements.txt

# Now copy the rest of the application files into the container
COPY ./blast_query /app/blast_query

# Set the command to run the FastAPI app
ENTRYPOINT ["uvicorn", "blast_query.api:app"]
