import datetime
import logging
from typing import Dict, Any

from pythonjsonlogger import jsonlogger

from infrastructure.settings import settings


class NoLabsFormatter(jsonlogger.JsonFormatter):
    def add_fields(self, log_record: Dict[str, Any], record: logging.LogRecord, message_dict: Dict[str, Any]) -> None:
        super().add_fields(log_record, record, message_dict)

        log_record['ts'] = (
                datetime.datetime.utcnow().isoformat(timespec='milliseconds') + 'Z'
        )
        log_record['level'] = settings.logging_level
        log_record['env'] = settings.environment