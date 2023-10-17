from src.server import settings
from src.server.api_handlers import AminoAcidLabApiMockHandler, DrugTargetApiMockHandler, AminoAcidLabApiHandler, \
    DrugTargetApiHandler


def api_handlers_factory():
    is_test = settings.is_test

    if is_test:
        return AminoAcidLabApiMockHandler(), DrugTargetApiMockHandler()

    return AminoAcidLabApiHandler(), DrugTargetApiHandler()