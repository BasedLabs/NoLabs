import os


def is_test():
    if 'reinvent_environment' not in os.environ:
        raise Exception('You must set reinvent_environment variable <test|prod> before running. ')
    return os.environ['reinvent_environment'] == 'test'
