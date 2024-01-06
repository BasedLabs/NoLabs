from typing import Any

from conformations.api_models import Integrators
from openmm.app import *
from openmm import *


def create_integrator(integrator: Integrators) -> Any:
    if integrator == integrator.langevin:
        return LangevinIntegrator
    if integrator == integrator.langevin:
        return LangevinMiddleIntegrator
    if integrator == integrator.nose_hoover:
        return NoseHooverIntegrator
    if integrator == integrator.brownian:
        return BrownianIntegrator
    if integrator == integrator.variable_verlet:
        return VariableVerletIntegrator