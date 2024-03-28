from dockstream.core.ligand.ligand_input_parser import LigandInputParser

from dockstream.utils.dockstream_exceptions import *

from dockstream.core.RDkit.RDkit_ligand_preparator import RDkitLigandPreparator
from dockstream.core.OpenEye.OpenEye_ligand_preparator import OpenEyeLigandPreparator
from dockstream.core.Corina.Corina_ligand_preparator import CorinaLigandPreparator
from dockstream.core.Schrodinger.Ligprep_ligand_preparator import LigprepLigandPreparator
from dockstream.core.OpenEyeHybrid.Omega_ligand_preparator import OmegaLigandPreparator

from dockstream.utils.enums.docking_enum import DockingConfigurationEnum
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.logging_enums import LoggingConfigEnum


def embed_ligands(smiles, pool_number, pool, logger, ligand_number_start=0):
    # enums
    _LE = LoggingConfigEnum()
    _LP = LigandPreparationEnum()
    _DE = DockingConfigurationEnum()

    # 1) load and parse the input, whether from command-line or from the configuration
    # note, that if "args.smiles" is not None, those smiles will be use irrespective of the configuration
    lig_inp_parser = LigandInputParser(smiles=smiles,
                                       ligand_number_start=ligand_number_start,
                                       **pool)
    list_ligands = lig_inp_parser.get_ligands()
    if len(list_ligands) == 0:
        raise LigandPreparationFailed("No smiles found in input.")
    logger.log(f"Loaded {len(list_ligands)} molecules.", _LE.DEBUG)

    # 2) do the embedding
    if pool[_LP.TYPE] == _LP.TYPE_RDKIT:
        prep = RDkitLigandPreparator(ligands=list_ligands, pool_number=pool_number, **pool)
    elif pool[_LP.TYPE] == _LP.TYPE_OPENEYE:
        prep = OpenEyeLigandPreparator(ligands=list_ligands, pool_number=pool_number, **pool)
    elif pool[_LP.TYPE] == _LP.TYPE_CORINA:
        prep = CorinaLigandPreparator(ligands=list_ligands, pool_number=pool_number, **pool)
    elif pool[_LP.TYPE] == _LP.TYPE_LIGPREP:
        prep = LigprepLigandPreparator(ligands=list_ligands, pool_number=pool_number, **pool)
    elif pool[_LP.TYPE] == _LP.TYPE_OMEGA:
        prep = OmegaLigandPreparator(ligands=list_ligands, pool_number=pool_number, **pool)
    else:
        raise LigandPreparationFailed("Type of pool is unknown.")

    # generate 3D coordinates (embed), if not using SDF input
    if _LP.INPUT_TYPE not in pool[_LP.INPUT].keys() or pool[_LP.INPUT][_LP.INPUT_TYPE].upper() != _LP.INPUT_TYPE_SDF:
        prep.generate3Dcoordinates()
    else:
        logger.log(f"As input is SDF, coordinate generation is skipped.", _LE.INFO)

    # 3) (optional) do alignment
    if _LP.ALIGN in pool.keys():
        prep.align_ligands()

    # 4) (optional) write the molecules to the disk
    if _LP.OUTPUT in pool.keys():
        prep.write_ligands(path=pool[_LP.OUTPUT][_LP.OUTPUT_CONFORMERPATH],
                           format=pool[_LP.OUTPUT][_LP.OUTPUT_FORMAT])

    # 5) save the ligands in the respective pool
    #    note, that a "pool" represents an embedded collection of molecules
    return prep
