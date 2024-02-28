import logging
from pythonjsonlogger import jsonlogger

from biobuddy.api_models import *

_logger = logging.getLogger()

_logger.setLevel(logging.DEBUG)

_logHandler = logging.StreamHandler()
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    @staticmethod
    def send_message_request(request: SendMessageToBioBuddyRequest):
        d = request.as_log_dict()
        _logger.info('Send message to biobuddy request', extra=d)

    @staticmethod
    def folding_response(response: SendMessageToBioBuddyResponse):
        d = response.as_log_dict()
        _logger.info('Send message to biobuddy response', extra=d)

    @staticmethod
    def exception():
        _logger.exception('Exception occurred in Bio buddy')

logger = Log()