from nolabs.application.use_cases.test1.workflow import Test1Component
from nolabs.application.use_cases.test2.workflow import Test2Component
from nolabs.application.workflow import Workflow

Workflow.set_component_types([Test1Component, Test2Component])