import sys
from pathlib import Path

import yaml


def parse_compose_service_build_context(microservice_name):
    compose = yaml.safe_load(Path("docker-compose.yaml").read_text())
    build_context = compose["services"][microservice_name]["build"]["context"]
    print(build_context)


if __name__ == "__main__":
    microservice_name = sys.argv[1]
    parse_compose_service_build_context(microservice_name)
