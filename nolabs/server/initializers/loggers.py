import logging
from nolabs.server.initializers.initialize_app import init

app, socketio = init()


class SocketIOEmitFilter(logging.Filter):
    def filter(self, record):
        if record.levelno == logging.INFO:
            # Emit the message using socketio
            socketio.emit('get-logs', {'response': record.msg})
        return True


class ConformationsSocketIOEmitFilter(logging.Filter):
    def filter(self, record):
        if record.levelno == logging.ERROR:
            # Emit the message using socketio
            socketio.emit('conformations-errors', {'response': record.msg})
        return True


logger = logging.getLogger('my_logger')
logger.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.addFilter(SocketIOEmitFilter())
logger.addHandler(ch)

conformations_errors_logger = logging.getLogger('conformations_logger')
conformations_errors_logger.setLevel(logging.DEBUG)

ch_conformations = logging.StreamHandler()
ch_conformations.setLevel(logging.DEBUG)
ch_conformations.addFilter(ConformationsSocketIOEmitFilter())
conformations_errors_logger.addHandler(ch_conformations)
