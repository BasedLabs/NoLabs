import sys
import yaml
from pathlib import Path


def parse_compose_service_image(microservice_name):
    compose = yaml.safe_load(Path('compose.yaml').read_text())
    return compose['services'][microservice_name]['image']


if __name__ == '__main__':
    microservice_name = sys.argv[1]
