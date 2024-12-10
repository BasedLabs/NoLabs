FROM python:3.12
WORKDIR /app
COPY ./requirements.txt /app/
RUN pip install --upgrade pip
RUN pip install --default-timeout=100 -r requirements.txt
COPY service /app