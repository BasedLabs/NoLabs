import dataclasses
from typing import List

import pydantic.dataclasses
from malevich.square import scheme
import localisation.api_models as nl


@scheme()
@dataclasses.dataclass
@pydantic.dataclasses.dataclass
class RunLocalisationPredictionRequest:
    amino_acid_sequence: str

    def map_to(self, model: 'RunLocalisationPredictionRequest') -> nl.RunLocalisationPredictionRequest:
        return nl.RunLocalisationPredictionRequest(
            amino_acid_sequence=model.amino_acid_sequence
        )


@scheme()
@dataclasses.dataclass
@pydantic.dataclasses.dataclass
class RunLocalisationPredictionResponse:
    cytosolic_proteins: float
    mitochondrial_proteins: float
    nuclear_proteins: float
    other_proteins: float
    extracellular_secreted_proteins: float
    errors: List[str]

    @staticmethod
    def map_from(model: nl.RunLocalisationPredictionResponse) -> 'RunLocalisationPredictionResponse':
        return RunLocalisationPredictionResponse(
            cytosolic_proteins=model.cytosolic_proteins,
            mitochondrial_proteins=model.mitochondrial_proteins,
            nuclear_proteins=model.nuclear_proteins,
            other_proteins=model.other_proteins,
            extracellular_secreted_proteins=model.extracellular_secreted_proteins,
            errors=model.errors
        )