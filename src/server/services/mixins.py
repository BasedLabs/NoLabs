import uuid


class UUIDGenerator():
    def gen_uuid(self):
        return {'id': str(uuid.uuid4())}