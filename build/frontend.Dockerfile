FROM node:lts

COPY . /app
WORKDIR /app/frontend

RUN npm install