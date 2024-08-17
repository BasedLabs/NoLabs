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