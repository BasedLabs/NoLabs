__all__ = [
    'mongo_connect'
]


from mongoengine import connect


def mongo_connect(connection_string: str):
    connect(host=connection_string)