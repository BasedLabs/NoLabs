FROM apache/airflow:2.9.3 AS airflow-stage

FROM python:3.11-buster

COPY --from=airflow-stage /opt/airflow /opt/airflow
COPY --from=airflow-stage /home/airflow /home/airflow

ENV AIRFLOW_HOME=/opt/airflow
ENV AIRFLOW_USER_HOME_DIR=/home/airflow

RUN pip install "apache-airflow[celery,redis]==2.9.3" --constraint "https://raw.githubusercontent.com/apache/airflow/constraints-2.9.3/constraints-3.8.txt"
RUN pip install psycopg2-binary

RUN useradd -ms /bin/bash airflow
RUN chown -R airflow: $AIRFLOW_HOME
USER airflow

RUN apt-get update; apt-get install curl gpg -y; \
mkdir -p /etc/apt/keyrings; \
curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg; \
echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_20.x nodistro main" | tee /etc/apt/sources.list.d/nodesource.list; \
 apt-get update && apt-get install -y nodejs;

RUN pip install poetry

COPY . /app
WORKDIR /app/frontend
RUN npm install
RUN npm run build

ENV PYTHONBUFFERED=1
ENV COMPOSE=1

RUN poetry install --no-root
WORKDIR /app
USER root
RUN chmod +x build/start_server.sh

RUN pip install uvicorn

RUN ls -la
ENTRYPOINT ["sh", "./build/start_server.sh"]