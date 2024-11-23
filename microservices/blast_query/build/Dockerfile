FROM python:3.11
RUN apt-get update && apt-get install -y curl python3-pip
WORKDIR /app
COPY ./requirements.txt /app/
RUN pip install --upgrade pip
RUN pip install --default-timeout=100 -r requirements.txt
COPY service /app
