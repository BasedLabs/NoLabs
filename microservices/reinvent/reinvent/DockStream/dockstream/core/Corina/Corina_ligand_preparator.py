import os
import shutil
from copy import deepcopy
from typing import Tuple, Optional, List
from typing_extensions import Literal
import tempfile

from pydantic import BaseModel, PrivateAttr
from rdkit import Chem

from dockstream.core.RDkit.RDkit_ligand_preparator import RDkitLigandPreparator
from dockstream.core.ligand.ligand import get_next_enumeration_number_for_ligand, Ligand, reset_enumerations_for_ligands

from dockstream.core.ligand_preparator import LigandPreparator, _LE

from dockstream.utils.dockstream_exceptions import LigandPreparationFailed
from dockstream.utils.general_utils import gen_temp_file

from dockstream.utils.smiles import to_mol
from dockstream.utils.execute_external.Corina import CorinaExecutor
from dockstream.utils.enums.Corina_enums import CorinaLigandPreparationEnum, CorinaExecutablesEnum

_LP = CorinaLigandPreparationEnum()
_EE = CorinaExecutablesEnum()


class Parallelization(BaseModel):
    number_cores: int = 1
    max_compounds_per_subjob: Optional[int] = None


class CorinaLigandPreparatorParameters(BaseModel):
    prefix_execution: Optional[str] = None
    binary_location: Optional[str] = None
    parallelization: Optional[Parallelization] = Parallelization()
    enumerate_stereo: Optional[bool] = False
    d_options: Optional[List[str]] = ["wh", "stergen", "preserve",
                                      "noflapn", "ori", "ampax",
                                      "names", "rc", "mc=1"]


class CorinaLigandPreparator(LigandPreparator, BaseModel):
    """Class that acts as an interface to the "Corina" executable to prepare ligands."""

    type: Literal["Corina"] = "Corina"
    parameters: CorinaLigandPreparatorParameters = CorinaLigandPreparatorParameters()

    class Config:
        underscore_attrs_are_private = True

    _Corina_executor: CorinaExecutor = PrivateAttr()

    def __init__(self, **data):
        super().__init__(**data)

        self._Corina_executor = CorinaExecutor(prefix_execution=self.parameters.prefix_execution,
                                               binary_location=self.parameters.binary_location)
        if not self._Corina_executor.is_available():
            raise LigandPreparationFailed("Cannot initialize Corina backend - abort.")
        self._logger.log(f"Checked Corina backend availability (prefix_execution={self.parameters.prefix_execution}, binary_location={self.parameters.binary_location}).",
                         _LE.DEBUG)

    def _load_references(self):
        references = []
        ref_format = self.align.reference_format.upper()
        for path in self.align.reference_paths:
            if ref_format == _LP.ALIGN_REFERENCE_FORMAT_PDB:
                ref_mol = Chem.MolFromPDBFile(path, sanitize=True)
                ref_mol.SetProp("_Name", os.path.basename(path))
                references.append(ref_mol)
            elif ref_format == _LP.ALIGN_REFERENCE_FORMAT_SDF:
                mol_supplier = Chem.SDMolSupplier(path)
                for mol in mol_supplier:
                    references.append(mol)
            else:
                raise IOError("Specified format not supported.")
        if len(references) == 0:
            raise LigandPreparationFailed("No reference molecules could be loaded with path(s) specified.")
        self._references = references
        self._logger.log(f"Stored {len(references)} reference molecules.", _LE.DEBUG)

    def _get_RDkit_aligner(self, conf, ligands):
        return RDkitLigandPreparator(ligands=ligands, **conf)

    def _get_d_parameters(self) -> str:
        if isinstance(self.parameters.d_options, str):
            self.parameters.d_options = [self.parameters.d_options]
        if not isinstance(self.parameters.d_options, list) or len(self.parameters.d_options) == 0:
            err_msg = f"If specified, parameter {_LP.D_OPTIONS} must be a list of strings."
            self._logger.log(err_msg, _LE.ERROR)
            raise ValueError(err_msg)
        d_parameters = ','.join(self.parameters.d_options)
        return d_parameters

    def _parse_molecules(self, tmp_sdf_path: str) -> List[Ligand]:
        mol_supplier = Chem.SDMolSupplier(tmp_sdf_path, removeHs=False)
        expanded_ligands = []
        for mol in mol_supplier:
            # Corina has a strange way of naming the conformers (e.g. "0:0_i001_c001" and "0:0_i002_c001" are
            # conformers of the same ligand ordered by internal strain energy; also, the option "mc=1" does not
            # reduce really to if multiple stereo-isomers are given
            if mol is not None and mol.HasProp("_Name"):
                name_parts = mol.GetProp("_Name").split('_')

                # check, that only one conformation per enumeration is taken forward
                if name_parts[2] != "c001":
                    continue
                for lig in self.ligands:
                    if name_parts[0] == lig.get_identifier():
                        # check, if it is the first (energy minimized) one and add it in case
                        # stereo-enumeration is disabled
                        if self.parameters.enumerate_stereo or name_parts[1] == "i001":
                            expanded_ligands.append(Ligand(smile=Chem.MolToSmiles(mol, isomericSmiles=True),
                                                           original_smile=lig.get_original_smile(),
                                                           ligand_number=lig.get_ligand_number(),
                                                           enumeration=lig.get_enumeration(),
                                                           molecule=mol,
                                                           mol_type=_LP.TYPE_CORINA,
                                                           name=lig.get_name()))
            else:
                self._logger.log(
                    "Skipped molecule when loading as _Name property could not be found - typically, this indicates that Corina could not embed the molecule.",
                    _LE.DEBUG)
        return expanded_ligands

    def _smiles_to_molecules(self, ligands: List[Ligand]) -> List[Ligand]:
        for lig in ligands:
            mol = to_mol(lig.get_smile())
            lig.set_molecule(mol)
            lig.set_mol_type(_LP.TYPE_CORINA)
        return ligands

    def generate3Dcoordinates(self):
        for lig in self.ligands:
            lig.set_molecule(None)
            lig.set_mol_type(None)
        ligand_list = self._smiles_to_molecules(deepcopy(self.ligands))

        # 1) generate temporary folder and files
        tmp_output_dir_path = tempfile.mkdtemp()
        tmp_smiles_path = gen_temp_file(suffix=".smi", dir=tmp_output_dir_path)
        tmp_molecules_path = gen_temp_file(suffix=".sdf", dir=tmp_output_dir_path)

        # 2) save the SMILES
        with open(tmp_smiles_path, 'w') as f:
            for lig in ligand_list:
                f.write(lig.get_smile() + " " + lig.get_identifier() + "\n")

        # 3) get "-d" parameters (either default or user specified)
        d_parameters = self._get_d_parameters()

        # 4) run "Corina" backend
        result = self._Corina_executor.execute(command=_EE.CORINA,
                                               arguments=[_EE.CORINA_D, d_parameters,
                                                          _EE.CORINA_T, _EE.CORINA_T_DISABLED,
                                                          _EE.CORINA_I, _EE.CORINA_T_SMILES, tmp_smiles_path,
                                                          _EE.CORINA_O, _EE.CORINA_T_SDF, tmp_molecules_path],
                                               check=False)
        self._logger.log(f"Executed Corina backend (output file: {tmp_molecules_path}).", _LE.DEBUG)

        # 5) load and store the conformers; name it sequentially
        #    note, that some backends require all H-coordinates (such as Glide) - so keep them!
        expanded_ligands = self._parse_molecules(tmp_molecules_path)

        # 6) merge newly embedded ligands with the old list
        merged_list = []
        for lig_old in self.ligands:
            # make a list with all the new enumerations for a given "old" ligand
            lig_enums_list = [lig_enum for lig_enum in expanded_ligands if lig_enum.get_identifier() == lig_old.get_identifier()]
            if len(lig_enums_list) == 0:
                # embedding failed completely, keep the old ligand (with "molecule" set to "None")
                merged_list.append(lig_old)
            else:
                # embedding succeeded, so replace the original ligand with the one (or more) embedded enumerations
                for lig_emb in lig_enums_list:
                    merged_list.append(lig_emb)
        reset_enumerations_for_ligands(merged_list)
        self.ligands = merged_list

        not_embedded = len([True for lig in self.ligands if lig.get_molecule() is None])
        if not_embedded > 0:
            self._logger.log(f"Corina might have had issues embedding all {len(self.ligands)} ligands, {not_embedded} were not obtained.",
                             _LE.WARNING)
            for lig in self.ligands:
                if lig.get_molecule() is None:
                    self._logger.log(f"It appears, Corina could not embed ligand {lig.get_identifier()} (smile: {lig.get_smile()}).",
                                     _LE.DEBUG)

        # 7) remove temporary files
        shutil.rmtree(tmp_output_dir_path)
        self._logger.log(f"In total, {len([True for lig in self.ligands if lig.get_molecule() is not None])} ligands (including enumerations) embedded (Corina backend).", _LE.DEBUG)

    def align_ligands(self):
        self.ligands = self._align_ligands_with_RDkit_preparator(self.ligands)
