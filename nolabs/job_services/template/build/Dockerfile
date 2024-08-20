ARG AIRFLOW_VERSION=2.9.3
ARG AIRFLOW_ENABLED=False

FROM apache/airflow:${AIRFLOW_VERSION} AS airflow-stage

FROM python:3.8

COPY --from=airflow-stage /opt/airflow /opt/airflow
COPY --from=airflow-stage /home/airflow /home/airflow

ENV AIRFLOW_HOME=/opt/airflow
ENV AIRFLOW_UID="1000"
ENV AIRFLOW_USER_HOME_DIR=/home/airflow

RUN apt-get update && apt-get install -y curl python3-pip

WORKDIR /app

COPY ./requirements.txt /app/

RUN python3.8 -m pip install --upgrade pip
RUN python3.8 -m pip install --default-timeout=100 -r requirements.txt

COPY ./esmfold_light /app/esmfold_light

ENTRYPOINT ["uvicorn", "esmfold_light.api:app"]
