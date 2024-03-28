from copy import deepcopy

from typing import List, Optional, Dict, Union
from pydantic import BaseModel, PrivateAttr
from rdkit import Chem, RDLogger
from dockstream.core.RDkit.RDkit_stereo_enumerator import RDKitStereoEnumerator

from dockstream.core.ligand.ligand import find_ligand

from dockstream.core.TautEnum.taut_enum_smile_preparation import TautEnumSmilePreparator
from dockstream.core.factories.transformator_factory import TransformatorFactory

from dockstream.loggers.ligand_preparation_logger import LigandPreparationLogger
from dockstream.utils.dockstream_exceptions import LigandPreparationFailed

from dockstream.utils.enums.docking_enum import DockingConfigurationEnum
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.RDkit_enums import RDkitLigandPreparationEnum
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.utils.enums.transformations_enums import TransformationEnum
from dockstream.utils.enums.stereo_enumeration_enums import StereoEnumerationEnum
from dockstream.utils.files_paths import generate_folder_structure

_DE = DockingConfigurationEnum()
_LP = LigandPreparationEnum()
_RLP = RDkitLigandPreparationEnum()
_LE = LoggingConfigEnum()
_TE = TransformationEnum()
_SE = StereoEnumerationEnum()


class CSVInput(BaseModel):
    smiles: str
    names: Optional[str] = None


class TautEnumInput(BaseModel):
    prefix_execution: str = None
    binary_location: str = None
    enumerate_protonation: bool = False


class AlignInput(BaseModel):
    mode: str
    reference_paths: List[str]
    reference_format: str
    minimum_substructure_ratio: float = 0.2
    fail_action: str = "keep"
    complete_rings_only: bool = True
    tethering: bool = False


class TransformationInput(BaseModel):
    type: str
    backend: str
    smirks: Optional[str]
    fail_action: str = "keep"


AnyStereoEnumerator = Union[RDKitStereoEnumerator]  # Add more when available.


class Input(BaseModel):
    type: Optional[str]
    input_path: Optional[str]
    tags: Optional[Dict[str, str]] = None
    delimiter: Optional[str] = ','
    initialization_mode: Optional[str] = _LP.INITIALIZATION_MODE_ORDER
    columns: Optional[CSVInput] = None
    use_taut_enum: Optional[TautEnumInput] = None
    stereo_enumeration: Optional[AnyStereoEnumerator] = None
    transformations: Optional[List[TransformationInput]] = None


class Output(BaseModel):
    format: str
    conformer_path: str


class LigandPreparator(BaseModel):
    """Base class implementing the interface for all docking preparation classes."""

    pool_id: str
    input: Input
    align: Optional[AlignInput] = None
    output: Optional[Output]
    ligands: Optional[List] = None

    _logger = PrivateAttr()
    _references: List = PrivateAttr(default=None)

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **data):
        super().__init__(**data)
        self._logger = LigandPreparationLogger()
        if self.ligands is not None and len(self.ligands) >= 1:
            self._initialize_ligands()

    def add_ligands(self, ligands):
        self.ligands = ligands
        self._initialize_ligands()

    def _initialize_ligands(self):
        # store ligands as list of "Ligands", generated either from smiles or molecules
        if not isinstance(self.ligands, list):
            self.ligands = [self.ligands]
        if len(self.ligands) == 0:
            raise LigandPreparationFailed("Specify at least one ligand (or a list).")

        # enumerate ligand smiles with tautomers / protomers, if specified
        if self.input.use_taut_enum is not None:
            self._taut_enum()

        # enumerate ligand smiles stereochemically, if specified
        if self.input.stereo_enumeration is not None:
            self._enumerate_stereoisomers()

        # apply transformations (e.g. SMIRKS), if specified
        if self.input.transformations is not None:
            self._apply_transformations()

        # treat the reference molecule(s) and store it internally as a list, if specified
        if self.align is not None:
            self._load_references()

    def _enumerate_stereoisomers(self):
        length_before = self.get_number_ligands()
        self.ligands = self.input.stereo_enumeration.enumerate(self.ligands)
        self._logger.log(f"Enumerated stereo-isomers (expanded {length_before} to {self.get_number_ligands()} enumerations).",
                         _LE.DEBUG)

    def _taut_enum(self):
        taut_enum = TautEnumSmilePreparator(enumerate_protonation=self.input.use_taut_enum.enumerate_protonation,
                                            original_enumeration=True,
                                            add_numbers_to_name=True,
                                            prefix_execution=self.input.use_taut_enum.prefix_execution,
                                            binary_location=self.input.use_taut_enum.binary_location)

        # taut_enum will return a list of "Ligand" objects, conditionally expanded by enumerated versions
        self.ligands = taut_enum.annotate_tautomers(ligands=self.ligands)
        self._logger.log("Executed taut_enum.", _LE.INFO)
        self._logger.log(f"Stored {len(self.ligands)} smiles from taut_enum output.", _LE.DEBUG)

    def _apply_transformations(self):
        number_transformations = 0
        list_transformators = TransformatorFactory(self.input.transformations).get_transformators()
        for transformator in list_transformators:
            self.ligands = transformator.transform(self.ligands)
        number_transformations += 1
        self._logger.log(f"Executed {number_transformations} transformation(s).", _LE.DEBUG)
        self._logger.log(f"After transformation stage, {len(self.ligands)} smiles were stored.", _LE.DEBUG)

    def set_references(self, references):
        # usually, references are loaded from files; but this function allows setting them as a list of molecules
        if references is not None:
            if not isinstance(references, list):
                references = [references]
        self._references = references

    def _load_references(self):
        raise NotImplementedError("This method is backend-specific and must be implemented by each individual child class.")

    def get_number_ligands(self):
        return len(self.ligands)

    def get_ligands(self):
        return self.ligands

    def get_number_references(self):
        if self._references is not None:
            return len(self._references)
        else:
            return None

    def get_references(self):
        return self._references

    def _get_RDkit_aligner(self, conf, ligands):
        raise NotImplementedError("This method is backend-specific and must be implemented by each individual child class.")

    def generate3Dcoordinates(self):
        raise NotImplementedError("This method is backend-specific and must be implemented by each individual child class.")

    def align_ligands(self):
        raise NotImplementedError("This method is backend-specific and must be implemented by each individual child class.")

    def _align_ligands_with_RDkit_preparator(self, ligands: list):
        if self.align.mode != _LP.ALIGN_MODE_INTERNAL:
            raise LigandPreparationFailed("Only internal alignment supported at the moment.")
        if self._references is None:
            raise LigandPreparationFailed("No reference molecule has been found.")

        # at this stage, "generate3Dcoordinates()" has been used to generate the conformers; use the internal alignment to a reference
        RDLogger.DisableLog("rdApp.*")
        rdkit_conf = {_LP.POOLID: "dummyPool",
                      _LP.INPUT: {},
                      _LP.TYPE: _LP.TYPE_RDKIT,
                      _LP.PARAMS: {
                          _RLP.EP_PARAMS_COORDGEN: {
                              _RLP.EP_PARAMS_COORDGEN_METHOD: _RLP.EP_PARAMS_COORDGEN_UFF,
                              _RLP.EP_PARAMS_COORDGEN_UFF_MAXITERS: 600
                          }
                      },
                      _LP.ALIGN: deepcopy(self.align.dict())}
        aligner = self._get_RDkit_aligner(conf=rdkit_conf,
                                          ligands=[deepcopy(lig) for lig in ligands])
        aligner.align_ligands()

        # overwrite molecules for those ligands, that could be aligned
        for aligned_lig in aligner.get_ligands():
            ligand = find_ligand(ligands=ligands,
                                 ligand_id=aligned_lig.get_ligand_number(),
                                 enumeration=aligned_lig.get_enumeration())
            if ligand is not None:
                ligand.set_molecule(aligned_lig.get_molecule())
        return ligands

    def write_ligands(self, path, format):
        format = format.upper()
        ligands_copy = [deepcopy(lig) for lig in self.ligands]

        # generate folder structure, if not available
        generate_folder_structure(filepath=path)

        # check and specify format of file
        # RDkit does not support the write-out of MOL2 files (apparently because of the format's inherent ambiguity)
        if format == _LP.OUTPUT_FORMAT_SDF:
            writer = Chem.SDWriter(path)
            for lig in ligands_copy:
                lig.add_tags_to_molecule()
                if lig.get_molecule() is not None:
                    mol = deepcopy(lig.get_molecule())
                    mol.SetProp("_Name", lig.get_identifier())
                    writer.write(mol)
            writer.close()
        elif format == _LP.OUTPUT_FORMAT_MAE:
            raise LigandPreparationFailed("Write-out as maestro file not yet implemented.")
        else:
            raise LigandPreparationFailed("Specified output format unknown.")
        self._logger.log(f"Wrote {len(self.ligands)} molecules to file {path} (format: {format}).", _LE.DEBUG)
