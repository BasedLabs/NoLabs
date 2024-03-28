from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum


class CorinaLigandPreparationEnum(LigandPreparationEnum):

    OUTPUT_FORMAT_MAE = "MAE"
    D_OPTIONS = "d_options"
    ENUMERATE_STEREO = "enumerate_stereo"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class CorinaOutputEnum:
    """This "Enum" serves to store all keywords that are used by the "corina" executable."""

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class CorinaExecutablesEnum:
    """This "Enum" serves to store all the executables (and parameters) as strings available in the "Corina" backend."""

    # executable "corina" + parameters
    # ---------
    CORINA = "corina"
    CORINA_I = "-i"
    CORINA_T_SMILES = "t=smiles"
    CORINA_O = "-o"
    CORINA_T_SDF = "t=sdf"
    CORINA_D = "-d"
    CORINA_HELP = "-h"
    CORINA_HELP_IDENTIFICATION_STRING = "3D-structure generator"     # if string found in "stderr" of result, "corina"
                                                                     # is available
    CORINA_T = "-t"                                                  # trace level (results in "corina.trc" file)
    CORINA_T_DISABLED = "n"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
