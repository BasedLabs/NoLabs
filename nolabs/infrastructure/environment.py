import argparse

__all__ = ['is_dev', ]

import os

if 'NOLABS_ENVIRONMENT' not in os.environ:
    raise Exception('You must specify NOLABS_ENVIRONMENT environment variable. Allowed values are "dev", "prod".')


def is_dev() -> bool:
    return os.environ.get('NOLABS_ENVIRONMENT') == 'dev'
