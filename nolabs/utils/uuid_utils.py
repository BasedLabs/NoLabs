import uuid


class UuidUtils:
    def generate_uuid(self) -> str:
        return str(uuid.uuid4())