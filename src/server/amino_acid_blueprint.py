from flask import Blueprint
from flask import request

from src.server.api_handlers import AminoAcidLabApiHandler


def resolve_api_endpoints(api_handler: AminoAcidLabApiHandler):
    amino_acid_bp = Blueprint('amino-acid', __name__)

    @amino_acid_bp.route('/inference', methods=['POST'])
    def inference_amino_acid():
        return api_handler.inference(request)

    return amino_acid_bp