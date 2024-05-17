class WorkflowException(Exception):
    def __init__(self, msg: str):
        super().__init__(msg)


class SchemaException(WorkflowException):
    def __init__(self, msg: str):
        super().__init__(msg)