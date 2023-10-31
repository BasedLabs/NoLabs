from flask_socketio import SocketIO

from flask import Flask
from flask_cors import CORS

app = None
socketio = None


def init():
    global app
    global socketio

    if not app:
        app = Flask(__name__)
        CORS(app)
        socketio = SocketIO(app, cors_allowed_origins='*')

    return app, socketio