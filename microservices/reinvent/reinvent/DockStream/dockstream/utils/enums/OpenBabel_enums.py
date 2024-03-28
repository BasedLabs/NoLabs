
class OpenBabelOutputEnum:

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class OpenBabelExecutablesEnum:

    # executable "obabel" + parameters
    # ---------
    OBABEL = "obabel"
    OBABEL_IDENTIFICATION_STRING = "-O<outfilename>"
    OBABLE_INPUTFORMAT_PDBQT = "-ipdbqt"                             # sets the input format to "PDBQT" (output of "AutoDock Vina")
    OBABEL_P = "-p"                                                  # sets the <pH> value (e.g. "-p 7.4") for protonation
                                                                     # note, that this overwrites "--addpolarh", which is thus not used
    OBABEL_O = "-O"                                                  # specifies the output path (directly pasted afterwards, e.g. "-Omypath.pdb")
    OBABEL_OUTPUT_FORMAT_PDBQT = "-opdbqt"                           # sets the output format to "PDBQT" (input for "AutoDock Vina")
    OBABEL_OUTPUT_FORMAT_SDF = "-osdf"                               # sets the output format to "SDF"
    OBABEL_X = "-x"                                                  # specifies generation options
    OBABEL_X_R = 'r'                                                 # one of the 'X' options ("-x"), which disables the tree construction of the receptor (makes it static), directly pasted together: e.g. "-xr"
    OBABEL_PARTIALCHARGE = "--partialcharge"                         # sets the partial charge generation method (execute "obabel -L charges" to see list of available methods)
    OBABEL_PARTIALCHARGE_GASTEIGER = "gasteiger"                     # one method to compute the partial charges, used as: "--partialcharge gasteiger"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
