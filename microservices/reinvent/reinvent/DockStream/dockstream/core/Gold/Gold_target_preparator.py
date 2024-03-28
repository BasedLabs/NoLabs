import os
import pickle
import ccdc
from ccdc.docking import Docker
from ccdc.io import MoleculeReader

from dockstream.core.target_preparator import TargetPreparator

from dockstream.utils.dockstream_exceptions import TargetPreparationFailed

from dockstream.utils.enums.Gold_enums import GoldTargetPreparationEnum, GoldTargetKeywordEnum
from dockstream.containers.target_preparation_container import TargetPreparationContainer


class GoldTargetPreparator(TargetPreparator):
    """Class that deals with all the target preparatory steps needed before docking using "GOLD" can commence."""

    def __init__(self, conf: TargetPreparationContainer, target, run_number=0):
        self._TP = GoldTargetPreparationEnum()
        self._TK = GoldTargetKeywordEnum()
        self._target_dict = {self._TK.VERSION: self._TK.CURRENT_VERSION}

        # invoke base class's constructor first
        super().__init__(conf=conf, run_number=run_number)

        # check, whether the backend run specified is an "GOLD" one
        if self._run_parameters[self._TP.RUNS_BACKEND] != self._TP.RUNS_BACKEND_GOLD:
            raise TargetPreparationFailed("Tried to make an GOLD preparation with different backend specification.")

        if isinstance(target, str):
            if os.path.isfile(target):
                _, file_extension = os.path.splitext(target)
                if file_extension == ".pdb":
                    self._target = Docker()
                    self._settings = self._target.settings
                    self._settings.add_protein_file(file_name=target)

                    # add information to dictionary
                    with open(target, 'r') as file:
                        self._target_dict[self._TK.TARGET_PDB] = [line for line in file]
                    self._target_dict[self._TK.TARGET_PDB_FILENAME] = os.path.basename(target)
                else:
                    raise TargetPreparationFailed("Specified input file must be in PDB format for GOLD.")
            else:
                raise TargetPreparationFailed("Input target file does not exist.")
        elif isinstance(target, ccdc.docking.Docker):
            raise NotImplementedError
        else:
            raise TargetPreparationFailed("Constructor only accepts and Protein.BindingSite object or a file path.")
        self._logger.log("Added target to GOLD settings.", self._TL.DEBUG)

    def specify_cavity(self):
        self._target_dict[self._TK.CAVITY_METHOD] = self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD]
        if self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD] == self._TP.CAVITY_METHOD_REFERENCE:
            # note, that "MoleculeReader" is able to discern many formats from the ending, including "mol2" and "pdb"
            ref_ligand_path = self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_PATH]
            ref_ligand = MoleculeReader(filename=ref_ligand_path)
            protein = self._settings.proteins[0]
            self._settings.binding_site = self._settings.BindingSiteFromLigand(protein, ref_ligand, distance=self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_DISTANCE])

            # add information to dictionary
            with open(ref_ligand_path, 'r') as file:
                self._target_dict[self._TK.REFERENCE_LIGAND] = [line for line in file]
            self._target_dict[self._TK.CAVITY_REFERENCE_DISTANCE] = self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_DISTANCE]
            self._target_dict[self._TK.REFERENCE_LIGAND_FILENAME] = os.path.basename(ref_ligand_path)
        elif self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD] == self._TP.CAVITY_METHOD_POINT:
            raise NotImplementedError
            # origin (x,x,x)
            # distance x
        else:
            raise TargetPreparationFailed("Specified cavity determination method not defined for GOLD.")
        self._logger.log(f"Generated GOLD Protein.BindingSite with method {self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD]}.", self._TL.DEBUG)

    def write_target(self, path):
        _, file_extension = os.path.splitext(path)
        if file_extension != ".pkl":
            raise TargetPreparationFailed("Receptor files must end on .pkl.")
        if self._TK.CAVITY_METHOD not in self._target_dict:
            self._logger.log("Need to have executed specify_cavity before writing out result - will attempt this now.", self._TL.WARNING)
            self.specify_cavity()
        with open(path, "wb") as f:
            pickle.dump(self._target_dict, f)
        self._logger.log(f"Wrote binding site to file {path}.", self._TL.DEBUG)
