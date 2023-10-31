import sys
import logging
from io import StringIO

from src.server.initializers.initialize_app import init

app, socketio = init()

class SocketIOEmitFilter(logging.Filter):
    def filter(self, record):
        if record.levelno == logging.INFO:
            # Emit the message using socketio
            socketio.emit('get-logs', {'response': record.msg})
        return True

logger = logging.getLogger('my_logger')
logger.setLevel(logging.DEBUG)

# Create console handler
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# Add the SocketIOEmit filter to the handler
ch.addFilter(SocketIOEmitFilter())

# Add the handler to the logger
logger.addHandler(ch)

class InterceptedStdout(StringIO):
    def __init__(self, original_stdout):
        super().__init__()
        self.original_stdout = original_stdout

    def write(self, s):
        # Write to the original stdout
        self.original_stdout.write(s)
        # Also save in the current object
        super().write(s)
        socketio.emit('get-logs', {'response': s})

    def getvalue(self):
        return super().getvalue()

class intercept_stdout:
    def __init__(self):
        self._original_stdout = sys.stdout
        self._intercepted = InterceptedStdout(self._original_stdout)

    def __enter__(self):
        sys.stdout = self._intercepted
        return self._intercepted

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stdout = self._original_stdout
