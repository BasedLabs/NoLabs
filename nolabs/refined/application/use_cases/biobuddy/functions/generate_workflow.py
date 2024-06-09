import json
from typing import Dict, Any

from nolabs.domain.experiment import ExperimentId
from nolabs.api_models.biobuddy import FunctionCall, FunctionCallReturnData
from nolabs.modules.biobuddy.functions.base_function import BiobuddyFunction, FunctionParameterDefinition
from nolabs.infrastructure.settings import Settings

components = []
class GenerateWorkflowFunction(BiobuddyFunction):
    def __init__(self, settings: Settings):
        parameters = [
            FunctionParameterDefinition(name="workflow",
                                        type="string",
                                        required=True,
                                        description=f"Generate a JSON which would be used to construct a workflow graph. Components: {', '.join(component.name for component in components)}",
                                        items_type="string")
        ]
        super().__init__("create_workflow", "Generate a workflow function. Creates a JSON from which a workflow is constructed", parameters)
        self._settings = settings
        self._components = components

    def execute(self, experiment_id: ExperimentId, arguments: Dict[str, Any]) -> FunctionCall:
        return_json_string = arguments[self.parameters[0].name]

        workflow = json.loads(return_json_string)

        return FunctionCall(function_name="create_workflow", parameters=[],
                            data=FunctionCallReturnData(files=workflow))

