__all__ = [
    'mongo_connect'
]

from mongoengine import connect, disconnect, get_db, ConnectionFailure

connection = None


def mongo_connect(connection_string: str):
    connection = connect(host=connection_string)
    return connection


def mongo_disconnect():
    disconnect()


def is_connected():
    try:
        # Try to get the database object for the given alias
        db = get_db('nolabs')
        return True
    except ConnectionFailure:
        return False