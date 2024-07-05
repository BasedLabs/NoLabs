__all__ = [
    'mongo_connect'
]

from mongoengine import connect, disconnect


connection = None


def mongo_connect(connection_string: str):
    connection = connect(host=connection_string)
    return connection


def mongo_disconnect():
    disconnect()
