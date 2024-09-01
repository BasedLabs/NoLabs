from nolabs.application.use_cases.proteins.workflow import ProteinsComponent
from nolabs.application.use_cases.folding.workflow import EsmfoldLightComponent
from nolabs.application.workflow import Workflow

Workflow.set_component_types([ProteinsComponent, EsmfoldLightComponent])
