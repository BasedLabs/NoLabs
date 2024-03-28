
class StereoEnumerationEnum:

    # stereo-enumeration
    # ---------
    STEREO_ENUM = "stereo_enumeration"
    STEREO_ENUM_BACKEND = "stereo_backend"
    STEREO_ENUM_BACKEND_RDKIT = "RDkit"
    STEREO_ENUM_PARAMETERS = "parameters"

    # RDkit
    # ---------
    RDKIT_TRY_EMBEDDING = "try_embedding"
    RDKIT_UNIQUE = "unique"
    RDKIT_MAX_ISOMERS = "max_isomers"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
