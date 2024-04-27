import sys
import yaml
from pathlib import Path


def parse_compose_service_build_context(microservice_name):
    compose = yaml.safe_load(Path('compose.yaml').read_text())
    print(compose['services'][microservice_name]['build']['context'])


if __name__ == '__main__':
    microservice_name = sys.argv[1]
