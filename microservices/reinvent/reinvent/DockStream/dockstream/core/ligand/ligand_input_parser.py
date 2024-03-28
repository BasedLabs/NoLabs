import os
import pandas as pd
from typing import Optional, Any
from pydantic import BaseModel, PrivateAttr

from rdkit import Chem

from dockstream.core.ligand_preparator import Input
from dockstream.loggers.ligand_preparation_logger import LigandPreparationLogger
from dockstream.utils.dockstream_exceptions import LigandPreparationFailed
from dockstream.core.ligand.ligand import Ligand

from dockstream.utils.smiles import standardize_smiles

from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.utils.enums.docking_enum import DockingConfigurationEnum
from dockstream.utils.smiles import to_smiles

_DE = DockingConfigurationEnum()
_LP = LigandPreparationEnum()
_LE = LoggingConfigEnum()


class LigandInputParser(BaseModel):
    """This class is able to parse various input specifications and produces a list of Ligand objects."""

    smiles: Optional[Any]
    ligand_number_start: int = 0
    input: Input

    _logger = PrivateAttr()

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **data):
        super().__init__(**data)
        self._logger = LigandPreparationLogger()

        # extract parts from the configuration for convenience and try to infer input type, if not explicitly stated
        if self.input.type is None:
            self._logger.log(
                "Input type has not been explicitly specified, will attempt to infer it - this is not recommended.",
                _LE.WARNING)
            self.input.type = self._infer_input_type()
        self.input.type = self.input.type.upper()

        #self._do_standardize_smiles = nested_get(self._pool_parameters, [_LP.INPUT,
        #                                                                 _LP.INPUT_STANDARDIZE_SMILES],
        #                                         default=False)

    def get_ligands(self) -> list:
        if self.input.type == _LP.INPUT_TYPE_CONSOLE:
            return self._ligands_from_console()
        elif self.input.type == _LP.INPUT_TYPE_LIST:
            return self._ligands_from_smiles_list(self.smiles)
        elif self.input.type == _LP.INPUT_TYPE_SMI:
            return self._ligands_from_smi_file()
        elif self.input.type == _LP.INPUT_TYPE_CSV:
            return self._ligands_from_csv_file()
        elif self.input.type == _LP.INPUT_TYPE_SDF:
            return self._ligands_from_sdf_file()
        else:
            raise LigandPreparationFailed(f"Input file type {self.input.type} is not supported.", _LE.ERROR)

    def _ligands_from_console(self) -> list:
        ligand_smiles = self.smiles.split(';')
        return self._ligands_from_smiles_list(ligand_smiles)

    def _ligands_from_smi_file(self) -> list:
        if self.input.input_path is None:
            self._logger.log("When using SMI input, an input path has to be specified.", _LE.ERROR)
        with open(self.input.input_path) as f_input:
            ligand_smiles = f_input.readlines()
        ligand_smiles = [x.strip() for x in ligand_smiles]
        return self._ligands_from_smiles_list(ligand_smiles)

    def _ligands_from_smiles_list(self, smiles: list) -> list:
        #if self._do_standardize_smiles:
        #    smiles = self._standardize_smiles(smiles)
        return_list = []
        for number_smile, smile in enumerate(smiles):
            return_list.append(Ligand(smile=smile,
                                      original_smile=smile,
                                      ligand_number=number_smile + self.ligand_number_start,
                                      enumeration=0,
                                      molecule=None,
                                      mol_type=None,
                                      name=None))
        return return_list

    def _ligands_from_csv_file(self) -> list:
        if self.input.input_path is None:
            self._logger.log("When using CSV input, an input path has to be specified.", _LE.ERROR)
        if self.input.columns.smiles is None:
            self._logger.log("When using CSV input, a smiles column has to be specified.", _LE.ERROR)

        # load data and check
        data = pd.read_csv(self.input.input_path,
                           delimiter=self.input.delimiter)
        if self.input.columns.smiles not in list(data.columns):
            raise LigandPreparationFailed(f"Could not find column {self.input.columns.smiles} in input file {self.input.input_path} with columns {list(data.columns)}.")
        names_ligands = None
        if self.input.columns.names is not None and self.input.columns.names in list(data.columns):
            names_ligands = [str(x) for x in data[self.input.columns.names].tolist()]

        # generate ligands
        ligands = self._ligands_from_smiles_list([str(x) for x in data[self.input.columns.smiles].tolist()])
        if names_ligands is not None and len(ligands) == len(names_ligands):
            for ligand, name in zip(ligands, names_ligands):
                ligand.set_name(name)
        return ligands

    def _ligands_from_sdf_file(self) -> list:
        if self.input.input_path is None:
            self._logger.log("When using SDF input, an input path has to be specified.", _LE.ERROR)
        lig_container = []
        mol_supplier = Chem.SDMolSupplier(self.input.input_path, removeHs=False)
        for mol_id, mol in enumerate(mol_supplier):
            name = None
            if self.input.tags is not None and _LP.INPUT_SDF_TAGNAME_NAMES in self.input.tags.keys():
                name_tag = self.input.tags[_LP.INPUT_SDF_TAGNAME_NAMES]
                if mol.HasProp(name_tag):
                    name = str(mol.GetProp(name_tag))
                else:
                    self._logger.log(f"Molecule number {mol_id} in input SDF file does not have name tag {name_tag} - will set to None.",
                                     _LE.DEBUG)
            if self.input.initialization_mode == _LP.INITIALIZATION_MODE_ORDER:
                lig_container.append(Ligand(smile=to_smiles(mol),
                                            original_smile=to_smiles(mol),
                                            ligand_number=mol_id,
                                            molecule=mol,
                                            mol_type=_LP.TYPE_RDKIT,
                                            name=name))
            elif self.input.initialization_mode == _LP.INITIALIZATION_MODE_AZDOCK:
                # TODO: fix / handle case where docked poses (with X:X:X) are fed in
                parts = str(mol.GetProp("_Name")).split(':')
                lig_container.append(Ligand(smile=to_smiles(mol),
                                            original_smile=to_smiles(mol),
                                            ligand_number=int(parts[0]),
                                            enumeration=int(parts[1]),
                                            molecule=mol,
                                            mol_type=_LP.TYPE_RDKIT,
                                            name=name))
            else:
                raise ValueError(f"Initialization mode {self.input.initialization_mode} is not supported.")
        return lig_container

    def _standardize_smiles(self, smiles: list) -> list:
        # TODO: think about removing this altogether
        ligand_smiles = standardize_smiles(smiles, min_heavy_atoms=2, max_heavy_atoms=500,
                                           element_list=None, remove_long_side_chains=False,
                                           neutralise_charges=False)
        self._logger.log("Ligand smiles have been standardized.", _LE.DEBUG)
        return ligand_smiles

    def _infer_input_type(self) -> str:
        if self.smiles is not None or self.input.input_path is None:
            if isinstance(self.smiles, list):
                self._logger.log(f"Inferred pool type {_LP.INPUT_TYPE_LIST}.", _LE.WARNING)
                return _LP.INPUT_TYPE_LIST
            else:
                self._logger.log(f"Inferred pool type {_LP.INPUT_TYPE_CONSOLE}.", _LE.WARNING)
                return _LP.INPUT_TYPE_CONSOLE
        else:
            _, ext = os.path.splitext(self.input.input_path)
            ext = ext.lstrip('.').upper()
            if ext == _LP.INPUT_TYPE_SMI:
                self._logger.log(f"Inferred pool type {_LP.INPUT_TYPE_SMI}.", _LE.WARNING)
                return _LP.INPUT_TYPE_SMI
            elif ext == _LP.INPUT_TYPE_CSV:
                self._logger.log(f"Inferred pool type {_LP.INPUT_TYPE_CSV}.", _LE.WARNING)
                return _LP.INPUT_TYPE_CSV
            elif ext == _LP.INPUT_TYPE_SDF:
                self._logger.log(f"Inferred pool type {_LP.INPUT_TYPE_SDF}.", _LE.WARNING)
                return _LP.INPUT_TYPE_SDF
            else:
                raise LigandPreparationFailed("Could not make educated guess on input type - abort.")
