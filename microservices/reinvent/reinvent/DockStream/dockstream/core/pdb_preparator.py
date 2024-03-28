import warnings
from pdbfixer import PDBFixer
from openmm.vec3 import Vec3
from openmm.app import PDBFile

from dockstream.loggers.target_preparation_logger import TargetPreparationLogger
from dockstream.utils.dockstream_exceptions import TargetPreparationFailed

from dockstream.containers.target_preparation_container import TargetPreparationContainer
from dockstream.utils.enums.target_preparation_enum import TargetPreparationEnum
from dockstream.utils.enums.logging_enums import LoggingConfigEnum


class PDBPreparator:
    """Wrapper class for "PDBFixer" functionality."""

    def __init__(self, conf: TargetPreparationContainer):
        self._TE = TargetPreparationEnum()
        self._TL = LoggingConfigEnum()
        self._logger = TargetPreparationLogger()
        self._config = conf

    def fix_pdb(self, input_pdb_file, output_pdb_file):
        """Function that loads a PDB file and writes a fixed version to another PDB file."""

        if self._TE.FIX not in self._config[self._TE.TARGETPREP].keys():
            raise TargetPreparationFailed("Cannot fix target, if respective configuration block is missing.")

        paras = self._config[self._TE.TARGETPREP][self._TE.FIX]

        # load the target from a PDB file (generated temporarily in the child classes)
        fixer = PDBFixer(filename=input_pdb_file)
        self._logger.log(f"Initialized PDBFixer.", self._TL.DEBUG)

        # perform fixes specified in the configuration
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        if paras[self._TE.FIX_STANDARDIZE]:
            fixer.replaceNonstandardResidues()
            self._logger.log(f"Replaced {len(fixer.nonstandardResidues)} non-standard residues.", self._TL.DEBUG)
        if paras[self._TE.FIX_REMOVEHETEROGENS]:
            fixer.removeHeterogens(keepWater=True)
            self._logger.log("Removed heterogens.", self._TL.DEBUG)
        if paras[self._TE.FIX_MISSINGHEAVYATOMS]:
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            self._logger.log(f"Added {len(fixer.missingAtoms)} missing atoms.", self._TL.DEBUG)
        if paras[self._TE.FIX_MISSINGHYDROGENS]:
            fixer.addMissingHydrogens(pH=7.0)
            self._logger.log("Added missing hydrogens.", self._TL.DEBUG)
        if paras[self._TE.FIX_ADDWATERBOX]:
            # one could use the crystallographic unit cell, but that might be missing from the HEADER so go for a cubic
            # cell with some distance instead
            maxSize = max(max((pos[i] for pos in fixer.positions)) -
                          min((pos[i] for pos in fixer.positions)) for i in range(3))
            boxSize = maxSize * Vec3(1, 1, 1)
            fixer.addSolvent(boxSize)
            self._logger.log("Added water box.", self._TL.DEBUG)

        # write out the fixed target as another PDB file
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb_file, 'w'))
