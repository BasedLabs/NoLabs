from typing import List

from pydantic import dataclasses as pcdataclass


@pcdataclass.dataclass
class AminoAcidResponse:
    sequence: str
    name: str
    pdb_file_name: str
    pdb_file: str


@pcdataclass.dataclass
class RunFoldingResponse:
    experiment_id: str
    experiment_name: str
    amino_acids: List[AminoAcidResponse]


@pcdataclass.dataclass
class ExperimentFastaPropertyResponse:
    filename: str
    content: str


@pcdataclass.dataclass
class ExperimentPropertiesResponse:
    amino_acid_sequence: str | None
    fastas: List[ExperimentFastaPropertyResponse]


@pcdataclass.dataclass
class GetExperimentResponse:
    experiment_id: str
    experiment_name: str
    amino_acids: List[AminoAcidResponse]
    properties: ExperimentPropertiesResponse
