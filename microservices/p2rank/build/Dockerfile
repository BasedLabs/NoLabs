# Use Python 3.8 as the base image
FROM python:3

# Install default JDK and JRE
RUN apt-get update && apt-get install -y default-jdk default-jre

# Set a working directory
WORKDIR /app

# Copy the application files into the container
COPY ./p2rank /app/p2rank
COPY ./requirements.txt /app

# Install Python dependencies
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

# Download and extract p2rank
RUN curl -o p2rank.tar.gz -L "https://github.com/rdk/p2rank/releases/download/2.4.1/p2rank_2.4.1.tar.gz" \
    && mkdir -p /app/p2rank/p2rank_source \
    && tar -xzvf p2rank.tar.gz -C /app/p2rank/p2rank_source \
    && rm p2rank.tar.gz


# Set the command to run the FastAPI app
ENTRYPOINT ["uvicorn", "p2rank.api:app"]
