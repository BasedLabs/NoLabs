import sys
from pathlib import Path

import yaml


def parse_compose_service_dockerfile(microservice_name):
    compose = yaml.safe_load(Path("docker-compose.yaml").read_text())
    image = compose["services"][microservice_name]["build"]["dockerfile"]
    print(image)


if __name__ == "__main__":
    microservice_name = sys.argv[1]
    parse_compose_service_dockerfile(microservice_name)
