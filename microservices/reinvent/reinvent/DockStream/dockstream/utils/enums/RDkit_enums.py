from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum


class RDkitLigandPreparationEnum(LigandPreparationEnum):

    TAG_RDOCK_TETHERED_ATOMS = "TETHERED ATOMS"

    # 3D coordinates generation
    # ---------
    EP_PARAMS_COORDGEN = "coordinate_generation"
    EP_PARAMS_COORDGEN_METHOD = "method"
    EP_PARAMS_COORDGEN_UFF = "UFF"
    EP_PARAMS_COORDGEN_UFF_MAXITERS = "maximum_iterations"

    # tie molecules to a reference during docking
    ALIGN_TETHERING = "tethering"

    # keywords for molecule tags
    # ---------
    TAG_ALIGNED_ATOMS = "ALIGNED ATOMS"
    TAG_ALIGNED_RMSD = "ALIGNED RMSD"
    TAG_ALIGNED_REFERENCE = "ALIGNED REFERENCE"
    TAG_ALIGNED_ATOMS_RATIO = "ALIGNED ATOMS RATIO"

    PROTONATE = "protonate"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")