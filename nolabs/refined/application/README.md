## Application layer guide

One concern per one module inside `application`
Each context has:

- RunJobFeature - takes job_id as input and runs job, returns JobResponse

- SetupJobFeature - takes inputs (optional id, name, etc), experiment_id, and setups Job, returns JobResponse

- GetJobFeature - takes job_id as input and returns JobResponse

- GetJobStatusFeature - returns status of job

set_inputs on the job overrides previous results

Dependencies methods must be Feature name without 'Feature' suffix. GetJobFeature -> get_job
Controller methods names must be Feature name without 'Feature' suffix RunJobFeature -> run_job

Wrap features body in try catch with Exception -> NoLabsException wrapping

# TODO
- To set up domain models changes reflect jobs set_result
- To change set_input for jobs to `__init__`
- To add mappings between exceptions and return codes
- Logging
- Self code review