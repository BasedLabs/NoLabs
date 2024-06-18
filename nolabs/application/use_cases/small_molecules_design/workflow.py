import uuid
from typing import List, Type

from pydantic import BaseModel

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.application.use_cases.small_molecules_design.use_cases import GetJobFeature, \
    RunLearningStageJobFeature, GetJobSmilesFeature, RunSamplingStageJobFeature
from nolabs.domain.models.common import Protein, JobId, JobName, Ligand, LigandName, LigandId, \
    DesignedLigandScore, DrugLikenessScore
from nolabs.domain.models.small_molecules_design import SmallMoleculesDesignJob
from nolabs.infrastructure.di import InfrastructureDependencies
from nolabs.workflow.component import Component, JobValidationError


class SmallMoleculesDesignLearningInput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]


class SmallMoleculesDesignLearningOutputItem(BaseModel):
    protein: uuid.UUID
    ligands: List[uuid.UUID]


class SmallMoleculesDesignLearningOutput(BaseModel):
    protein_ligands_pairs: List[SmallMoleculesDesignLearningOutputItem]


class SmallMoleculesDesignLearningComponent(Component[SmallMoleculesDesignLearningInput, SmallMoleculesDesignLearningOutput]):
    name = 'Small molecules design'

    async def execute(self):
        if await self.jobs_setup_errors():
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Jobs are not valid')

        run_learning_job_feature = RunLearningStageJobFeature(api=InfrastructureDependencies.reinvent_microservice())
        get_smiles_feature = GetJobSmilesFeature(api=InfrastructureDependencies.reinvent_microservice())

        result = []

        job: SmallMoleculesDesignJob
        for job in self.jobs:
            await run_learning_job_feature.handle(job_id=job.id)
            smiles = await get_smiles_feature.handle(job_id=job.id)

            protein = job.protein
            ligand_ids = []

            for smi in smiles:
                id = uuid.uuid4()
                ligand = Ligand.create(
                    experiment=self.experiment,
                    name=LigandName(f'Generated {id}'),
                    id=LigandId(id),
                    smiles_content=smi.smiles
                )
                ligand.set_designed_ligand_score(DesignedLigandScore(smi.score))
                ligand.set_drug_likeness_score(DrugLikenessScore(smi.drug_likeness))
                ligand.add_binding(protein=protein)
                ligand.save()

                ligand_ids.append(ligand.id)

            result.append(SmallMoleculesDesignLearningOutputItem(
                protein=job.protein.id,
                ligands=ligand_ids
            ))

        self.output = SmallMoleculesDesignLearningOutput(
            protein_ligands_pairs=result
        )

    async def setup_jobs(self):
        self.jobs = []

        for protein_id in self.input.proteins_with_pdb:
            protein = Protein.objects.with_id(protein_id)

            job_id = JobId(uuid.uuid4())
            job_name = JobName(f'Small molecules design job for protein {protein.name}')

            job = SmallMoleculesDesignJob(
                id=job_id,
                name=job_name,
                experiment=self.experiment,
                protein=protein
            )

            job.save()

            self.jobs.append(job)

    async def jobs_setup_errors(self) -> List[JobValidationError]:
        validation_errors = []

        job: SmallMoleculesDesignJob
        for job in self.jobs:
            if not job.input_valid():
                validation_errors.append(
                    JobValidationError(
                        job_id=job.id,
                        msg=f'Job input is invalid. Setup job inputs manually'
                    )
                )

        return validation_errors

    @property
    def _input_parameter_type(self) -> Type[SmallMoleculesDesignLearningInput]:
        return SmallMoleculesDesignLearningInput

    @property
    def _output_parameter_type(self) -> Type[SmallMoleculesDesignLearningOutput]:
        return SmallMoleculesDesignLearningOutput

# ---

class SmallMoleculesDesignSamplingInput(BaseModel):
    small_molecules_design_learning_job_ids: List[uuid.UUID]


class SmallMoleculesDesignSamplingOutputItem(BaseModel):
    protein: uuid.UUID
    ligands: List[uuid.UUID]


class SmallMoleculesDesignSamplingOutput(BaseModel):
    items: List[SmallMoleculesDesignSamplingOutputItem]


class SmallMoleculesDesignSamplingComponent(Component[SmallMoleculesDesignSamplingInput, SmallMoleculesDesignSamplingOutput]):
    name = 'Small molecules design sampling'

    async def execute(self):
        if await self.jobs_setup_errors():
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Jobs are not valid')

        run_sampling_job_feature = RunSamplingStageJobFeature(api=InfrastructureDependencies.reinvent_microservice())
        get_smiles_feature = GetJobSmilesFeature(api=InfrastructureDependencies.reinvent_microservice())

        result = []

        job: SmallMoleculesDesignJob
        for job_id in self.input.small_molecules_design_learning_job_ids:
            job = SmallMoleculesDesignJob.objects.with_id(job_id)
            await run_sampling_job_feature.handle(job_id=job_id)
            smiles = await get_smiles_feature.handle(job_id=job_id)

            protein = job.protein
            ligand_ids = []

            for smi in smiles:
                id = uuid.uuid4()
                ligand = Ligand.create(
                    experiment=self.experiment,
                    name=LigandName(f'Generated {id}'),
                    id=LigandId(id),
                    smiles_content=smi.smiles
                )
                ligand.set_designed_ligand_score(DesignedLigandScore(smi.score))
                ligand.set_drug_likeness_score(DrugLikenessScore(smi.drug_likeness))
                ligand.add_binding(protein=protein)
                ligand.save()

                ligand_ids.append(ligand.id)

            result.append(SmallMoleculesDesignSamplingOutputItem(
                protein=job.protein.id,
                ligands=ligand_ids
            ))

        self.output = SmallMoleculesDesignSamplingOutput(
            items=result
        )

    async def setup_jobs(self):
        pass

    async def jobs_setup_errors(self) -> List[JobValidationError]:
        job: SmallMoleculesDesignJob

        api = InfrastructureDependencies.reinvent_microservice()

        errors = []

        for job in self.jobs:
            contig_result = api.get_config_api_reinvent_reinvent_config_id_get(config_id=str(job.id))
            if not contig_result.actual_instance.sampling_allowed:
                errors.append(
                    JobValidationError(
                        job_id=job.id,
                        msg='You must run learning stage to run sampling. This job is not ready'
                    )
                )

        return errors

    @property
    def _input_parameter_type(self) -> Type[SmallMoleculesDesignSamplingInput]:
        return SmallMoleculesDesignSamplingInput

    @property
    def _output_parameter_type(self) -> Type[SmallMoleculesDesignSamplingOutput]:
        return SmallMoleculesDesignSamplingOutput