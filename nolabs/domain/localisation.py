from typing import Any

import pydantic.dataclasses as pcdataclasses
from pydantic import model_validator


@pcdataclasses.dataclass
class LocalisationProbability:
    cytosolic_proteins: float
    mitochondrial_proteins: float
    nuclear_proteins: float
    other_proteins: float
    extracellular_secreted_proteins: float

    @model_validator(mode='after')
    @classmethod
    def check_inputs(cls, data: Any) -> Any:
        if not isinstance(data, LocalisationProbability):
            raise ValueError('Incorrect data type')

        probabilities_sum = 0.0
        probabilities_sum += data.mitochondrial_proteins
        probabilities_sum += data.cytosolic_proteins
        probabilities_sum += data.nuclear_proteins
        probabilities_sum += data.other_proteins
        probabilities_sum += data.extracellular_secreted_proteins

        if abs(probabilities_sum - 1.0) > 0.01:
            raise ValueError('Probabilities sum are not 1.0')

        return data
