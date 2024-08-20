ARG AIRFLOW_VERSION=2.9.3
ARG AIRFLOW_ENABLED=False

FROM apache/airflow:${AIRFLOW_VERSION} AS airflow-stage

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
FROM python:3.8

COPY --from=airflow-stage /opt/airflow /opt/airflow
COPY --from=airflow-stage /home/airflow /home/airflow

ENV AIRFLOW_HOME=/opt/airflow
ENV AIRFLOW_UID="1000"
ENV AIRFLOW_USER_HOME_DIR=/home/airflow

RUN apt-get update && apt-get install -y curl python3-pip

WORKDIR /app

COPY . /app

RUN ls -l

RUN pip install poetry
RUN poetry install

ENTRYPOINT ["uvicorn", "esmfold_light.api:app"]
