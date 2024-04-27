import os
import sys
import yaml
from pathlib import Path


def parse_compose_service_build_context(microservice_name):
    compose = yaml.safe_load(Path('compose.yaml').read_text())
    build_context = compose['services'][microservice_name]['build']['context']
    env_file = os.getenv('GITHUB_ENV')
    with open(env_file, 'a') as fh:
        fh.write(f'BUILD_CONTEXT={build_context}')


if __name__ == '__main__':
    microservice_name = sys.argv[1]
