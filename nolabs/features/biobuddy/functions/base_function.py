from typing import Dict, Any, List

from pydantic import dataclasses as pcdataclass

from nolabs.api_models.biobuddy import FunctionCall
from nolabs.domain.experiment import ExperimentId


@pcdataclass.dataclass
class FunctionParameterDefinition:
    name: str
    type: str
    required: bool
    description: str
    items_type: str = None  # only required for arrays


class BiobuddyFunction:
    def __init__(self, name: str, description: str, parameters: List[FunctionParameterDefinition]):
        self.name = name
        self.description = description
        self.parameters = parameters

    def execute(self, experiment_id: ExperimentId, arguments: Dict[str, Any]) -> FunctionCall:
        raise NotImplementedError("Subclasses must implement this method")

# Example derived class
