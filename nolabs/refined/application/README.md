1 domain context per module
Each context has

RunJobFeature - takes job_id as input and runs job, returns JobResponse

SetupJobFeature - takes inputs (optional id, name, etc), experiment_id, and setups Job, returning JobResponse

GetJobFeature - takes job_id as input and returns JobResponse

GetJobStatusFeature - returns status of job

set_inputs on the job overrides previous results

Dependencies methods must be Feature name without 'Feature' suffix
Controller methods names must be Feature name without 'Feature' suffix

Wrap features body in try catch with Exception -> NoLabsException wrapping

# TODO add experiment id to responses