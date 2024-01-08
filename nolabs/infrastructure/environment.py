import argparse

__all__ = ['is_dev', ]

_dev = 'dev'
_prod = 'prod'
_parser = argparse.ArgumentParser()
_parser.add_argument('-n',
                     '--environment',
                     metavar='=ENVIRONMENT',
                     type=str,
                     help='Environment',
                     choices=[_dev, _prod],
                     required=True)
_args = _parser.parse_args()


def is_dev() -> bool:
    args = _parser.parse_args()
    return args.environment == _dev
