FROM python:3.11-buster

RUN apt-get update; apt-get install curl gpg -y; \
mkdir -p /etc/apt/keyrings; \
curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg; \
echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_20.x nodistro main" | tee /etc/apt/sources.list.d/nodesource.list; \
 apt-get update && apt-get install -y nodejs;

RUN pip install poetry

COPY . /app
WORKDIR /app/frontend
RUN npm install

ENV PYTHONBUFFERED=1

RUN poetry install --no-root
WORKDIR /app
RUN chmod +x build/start_server.sh

ENTRYPOINT ["./build/start_server.sh"]