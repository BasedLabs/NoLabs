
class TransformationEnum:

    TRANSFORMATIONS = "transformations"
    TRANSFORMATION_TYPE = "type"
    TRANSFORMATION_TYPE_SMIRKS = "smirks"
    TRANSFORMATION_BACKEND = "backend"
    TRANSFORMATION_BACKEND_OPENEYE = "OpenEye"
    TRANSFORMATION_SMIRKS = "smirks"
    TRANSFORMATION_FAIL_ACTION = "fail_action"
    TRANSFORMATION_FAIL_ACTION_KEEP = "keep"
    TRANSFORMATION_FAIL_ACTION_DISCARD = "discard"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
