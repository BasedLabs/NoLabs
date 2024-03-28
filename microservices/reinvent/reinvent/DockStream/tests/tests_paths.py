import json
from dockstream.utils.files_paths import attach_root_path

# load the instantiated "config.json", holding the license key for OpenEye for example
try:
    with open(attach_root_path("dockstream/config/tests_config/config.json"), 'r') as f:
        MAIN_CONFIG = json.load(f)
except:
    MAIN_CONFIG = {}


class PATHS_1UYD:
    TARGET_APO_PDB = "tests/tests_data/1UYD/1UYD_apo.pdb"
    TARGET_APO_MOL2 = "tests/tests_data/1UYD/1UYD_apo.mol2"
    TARGET_LIGAND_PU8_PDB = "tests/tests_data/1UYD/1UYD_ligand_PU8.pdb"
    LIGAND_PU8_SDF = "tests/tests_data/1UYD/ligand_PU8.sdf"
    LIGAND_PU8_PDB = "tests/tests_data/1UYD/PU8.pdb"
    LIGAND_MISSING_PARTS_PDB = "tests/tests_data/1UYD/1UYD_apo_missing_parts.pdb"
    LIGANDS_SDF = "tests/tests_data/1UYD/ligands.sdf"
    LIGANDS_SDF_WITH_HYDROGENS = "tests/tests_data/1UYD/ligands_with_hydrogens.sdf"
    LIGANDS_SDF_ALIGNED = "tests/tests_data/1UYD/ligands_aligned.sdf"
    LIGANDS_SDF_ALIGNED_TETHERED = "tests/tests_data/1UYD/ligands_aligned_tethered.sdf"
    LIGANDS_SMILES_TXT = "tests/tests_data/1UYD/ligands_smiles.txt"
    LIGANDS_SMILES_SMI_TAUT_ENUM = "tests/tests_data/1UYD/ligands_smiles_taut_enum.smi"
    LIGANDS_SDF_CORINA = "tests/tests_data/Corina/ligands.sdf"
    LIGANDS_WITH_ENUMERATION_SDF = "tests/tests_data/1UYD/ligands_with_enumeration.sdf"
    LIGANDS_CSV = "tests/tests_data/1UYD/ligands.csv"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class PATH_RDOCK_EXAMPLES:
    CAVITY_PRM = "tests/tests_data/rDock/rbcavity_1UYD.prm"
    CAVITY_PRM_UPDATED = "tests/tests_data/rDock/rbcavity_1UYD_updated.prm"
    DOCKED_LIGANDS = "tests/tests_data/rDock/output_ligands_tethered.sd"
    CAVITY_AS = "tests/tests_data/rDock/rbcavity_1UYD_updated.as"
    CAVITY_GRID = "tests/tests_data/rDock/rbcavity_1UYD_updated_cav1.grd"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class PATH_OPENEYE_EXAMPLES:
    RECEPTOR = "tests/tests_data/OpenEye/1UYD_reference_receptor.oeb"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class PATH_OPENEYEHYBRID_EXAMPLES:
    RECEPTOR = "tests/tests_data/OpenEyeHybrid/1UYD_reference_receptor.oeb"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class PATH_SCHRODINGER_EXAMPLES:
    GRIDFILE = "tests/tests_data/Schrodinger/1UYD_grid.zip"
    CONSTRAINTS_GRIDFILE = "tests/tests_data/Schrodinger/1UYD_hbond_constraints.zip"
    LICADMIN_FILE = "tests/tests_data/Schrodinger/licadmin_output.txt"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class PATH_GOLD_EXAMPLES:
    TARGETFILE = "tests/tests_data/Gold/Gold_binding_site.pkl"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class PATH_AUTODOCKVINA_EXAMPLES:
    RECEPTOR = "tests/tests_data/AutoDockVina/1UYD_fixed.pdbqt"
    BACKEND_TESTS_FOLDER = "tests/junk/AutoDockVina_backend"
    TARGET_PREP_FOLDER = "tests/junk/AutoDockVina_target_prep"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
