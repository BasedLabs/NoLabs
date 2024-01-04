import os

__all__ = ['is_dev',]

environment_env_key = 'NOLABS_ENV'


def is_dev():
    return os.getenv(environment_env_key)
