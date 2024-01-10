# RunGromacsSimulationsRequest

RunGromacsSimulationsRequest(*args: Any, temperatureK: float = 273.15, frictionCoeff: float = 1.0, stepSize: float = 0.002, integrator: conformations.api_models.Integrators = <Integrators.langevin: 'LangevinIntegator'>, takeFrameEvery: int = 1000, totalFrames: int = 10000, top: str, gro: str)

## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**temperature_k** | **float** |  | [optional] [default to 273.15]
**friction_coeff** | **float** |  | [optional] [default to 1.0]
**step_size** | **float** |  | [optional] [default to 0.002]
**integrator** | [**Integrators**](Integrators.md) |  | [optional] 
**take_frame_every** | **int** |  | [optional] [default to 1000]
**total_frames** | **int** |  | [optional] [default to 10000]
**top** | **str** |  | 
**gro** | **str** |  | 

## Example

```python
from conformations_microservice.models.run_gromacs_simulations_request import RunGromacsSimulationsRequest

# TODO update the JSON string below
json = "{}"
# create an instance of RunGromacsSimulationsRequest from a JSON string
run_gromacs_simulations_request_instance = RunGromacsSimulationsRequest.from_json(json)
# print the JSON string representation of the object
print RunGromacsSimulationsRequest.to_json()

# convert the object into a dict
run_gromacs_simulations_request_dict = run_gromacs_simulations_request_instance.to_dict()
# create an instance of RunGromacsSimulationsRequest from a dict
run_gromacs_simulations_request_form_dict = run_gromacs_simulations_request.from_dict(run_gromacs_simulations_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


