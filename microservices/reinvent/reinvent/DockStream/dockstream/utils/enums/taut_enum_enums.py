
class TautEnumEnum:
    """This "Enum" serves to store all keywords that are used by the "taut_enum" executable."""

    TAUTENUM = "taut_enum"
    TAUTENUM_I = "-I"
    TAUTENUM_O = "-O"
    TAUTENUM_HELP = "--help"
    TAUTENUM_HELP_IDENTIFICATION_STRING = "Allowed Options"
    TAUTENUM_ENUM_PROTO = "--enumerate-protonation"
    TAUTENUM_ORI_ENUM = "--original-enumeration"
    TAUTENUM_ADD_NUMBERS = "--add-numbers-to-name"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
