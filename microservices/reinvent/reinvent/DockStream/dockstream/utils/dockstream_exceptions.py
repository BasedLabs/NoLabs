
class LigandPreparationFailed(Exception):
    pass


class ConfigParsingFailed(Exception):
    pass


class DockingRunFailed(Exception):
    pass


class TargetPreparationFailed(Exception):
    pass


class ResultParsingFailed(Exception):
    pass


class TransformationFailed(Exception):
    pass


def get_exception_message(e: Exception):
    if hasattr(e, "message"):
        return e.message
    else:
        return e
