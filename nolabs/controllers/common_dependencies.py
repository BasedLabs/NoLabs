from nolabs.infrastructure.settings import Settings
from nolabs.utils import DateTimeUtils, uuid_utils


def settings_dependency() -> Settings:
    return Settings()


def dt_utils_dependency() -> DateTimeUtils:
    return DateTimeUtils()


def generate_uuid_dependency() -> uuid_utils.UuidUtils:
    return uuid_utils.UuidUtils()
