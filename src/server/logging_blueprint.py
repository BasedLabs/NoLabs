from flask import Blueprint, request
from src.server.initializers import loggers
import src.server.settings as settings


def resolve_api_endpoints():
    logs_bp = Blueprint('logs', __name__)

    @logs_bp.route('/conformations/logs', methods=['POST'])
    def logs():
        j = request.get_json(force=True)
        if j['get-logs'] and settings.CONFORMATIONS_IN_PROCESS:
            loggers.logger.info(j['get-logs'])
        return 'Ok', 200

    @logs_bp.route('/conformations/errors', methods=['POST'])
    def errors():
        j = request.get_json(force=True)
        if j['conformations-errors'] and settings.CONFORMATIONS_IN_PROCESS:
            loggers.conformations_errors_logger.error(j['conformations-errors'])
        return 'Ok', 200

    return logs_bp
