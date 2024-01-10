# RunPdbSimulationsRequest

RunPdbSimulationsRequest(*args: Any, temperatureK: float = 273.15, frictionCoeff: float = 1.0, stepSize: float = 0.002, integrator: conformations.api_models.Integrators = <Integrators.langevin: 'LangevinIntegator'>, takeFrameEvery: int = 1000, totalFrames: int = 10000, pdbContent: str, forceField: str, waterForceField: conformations.api_models.WaterForceFields)

## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**temperature_k** | **float** |  | [optional] [default to 273.15]
**friction_coeff** | **float** |  | [optional] [default to 1.0]
**step_size** | **float** |  | [optional] [default to 0.002]
**integrator** | [**Integrators**](Integrators.md) |  | [optional] 
**take_frame_every** | **int** |  | [optional] [default to 1000]
**total_frames** | **int** |  | [optional] [default to 10000]
**pdb_content** | **str** |  | 
**force_field** | **str** |  | 
**water_force_field** | [**WaterForceFields**](WaterForceFields.md) |  | 

## Example

```python
from conformations_microservice.models.run_pdb_simulations_request import RunPdbSimulationsRequest

# TODO update the JSON string below
json = "{}"
# create an instance of RunPdbSimulationsRequest from a JSON string
run_pdb_simulations_request_instance = RunPdbSimulationsRequest.from_json(json)
# print the JSON string representation of the object
print RunPdbSimulationsRequest.to_json()

# convert the object into a dict
run_pdb_simulations_request_dict = run_pdb_simulations_request_instance.to_dict()
# create an instance of RunPdbSimulationsRequest from a dict
run_pdb_simulations_request_form_dict = run_pdb_simulations_request.from_dict(run_pdb_simulations_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


