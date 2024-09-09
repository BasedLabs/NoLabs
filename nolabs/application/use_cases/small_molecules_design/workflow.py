import os
import uuid
from typing import List, Type

import asyncio
from pydantic import BaseModel

from domain.exceptions import NoLabsException, ErrorCodes
from nolabs.application.use_cases.small_molecules_design.use_cases import RunLearningStageJobFeature, \
    GetJobSmilesFeature, GetJobStatusFeature
from nolabs.domain.models.common import Protein, JobId, JobName, Ligand, LigandName, DesignedLigandScore, DrugLikenessScore
from nolabs.domain.models.small_molecules_design import SmallMoleculesDesignJob
from nolabs.infrastructure.di import InfrastructureDependencies
from nolabs.application.workflow.component import Component, JobValidationError


class SmallMoleculesDesignLearningInput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]


class SmallMoleculesDesignLearningOutputItem(BaseModel):
    protein: uuid.UUID
    ligands: List[uuid.UUID]


class SmallMoleculesDesignLearningOutput(BaseModel):
    protein_ligands_pairs: List[SmallMoleculesDesignLearningOutputItem]


class SmallMoleculesDesignLearningComponent(
    Component[SmallMoleculesDesignLearningInput, SmallMoleculesDesignLearningOutput]):
    name = 'Small molecules design'

    async def execute(self):
        if await self.jobs_setup_errors():
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Jobs are not valid')

        run_learning_job_feature = RunLearningStageJobFeature(api=InfrastructureDependencies.reinvent_microservice())
        get_smiles_feature = GetJobSmilesFeature(api=InfrastructureDependencies.reinvent_microservice())
        is_job_running_feature = GetJobStatusFeature(api=InfrastructureDependencies.reinvent_microservice())

        result = []

        job: SmallMoleculesDesignJob
        for job in self.jobs:
            try:

                await run_learning_job_feature.handle(job_id=job.id)

                await asyncio.sleep(0.5)

                running_status = await is_job_running_feature.handle(job_id=job.id)

                while running_status.running:
                    await asyncio.sleep(0.5)
                    running_status = await is_job_running_feature.handle(job_id=job.id)

                smiles = await get_smiles_feature.handle(job_id=job.id)

                protein = job.protein

                def normalize_floats(f: float) -> str:
                    return str(f).replace('.', 'dot')

                for smi in smiles:
                    name = (str(job.protein.name) + '-binder-drug-likeness-'
                            + normalize_floats(smi.drug_likeness) +
                            '-score-' + normalize_floats(smi.score) + '-stage-' + smi.stage)
                    ligand = Ligand.create(
                        experiment=job.experiment,
                        name=LigandName(name),
                        smiles_content=smi.smiles
                    )
                    ligand.set_designed_ligand_score(DesignedLigandScore(smi.score))
                    ligand.set_drug_likeness_score(DrugLikenessScore(smi.drug_likeness))
                    ligand.save()
                    ligand.add_binding(protein=protein, confidence=smi.score)

                    job.ligands.append(ligand)

                    await job.save()

            finally:
                job.finished()
                await job.save()

            result.append(SmallMoleculesDesignLearningOutputItem(
                protein=job.protein.id,
                ligands=[l.id for l in job.ligands]
            ))

        self.output = SmallMoleculesDesignLearningOutput(
            protein_ligands_pairs=result
        )

    async def setup_jobs(self):
        api = InfrastructureDependencies.reinvent_microservice()

        for protein_id in self.input.proteins_with_pdb:
            protein = Protein.objects.with_id(protein_id)

            job_id = JobId(uuid.uuid4())

            job = SmallMoleculesDesignJob(
                id=job_id,
                name=JobName(f'Small molecules design for {protein.name}'),
                experiment=self.experiment
            )

            if not protein.pdb_content:
                raise NoLabsException(ErrorCodes.protein_pdb_is_empty, messages='Protein pdb content is undefined')

            job.set_protein(protein=protein)

            job.set_inputs(
                protein=protein,
                center_x=0,
                center_y=0,
                center_z=0,
                size_x=5.0,
                size_y=5.0,
                size_z=5.0,
                batch_size=50,
                minscore=0.4,
                epochs=128,
                throw=False
            )

            tmp_file_path = 'tmp.pdb'
            open(tmp_file_path, 'wb').write(job.protein.pdb_content)

            api.save_params_api_reinvent_config_id_params_post(
                config_id=str(job.iid.value),
                name=job.name.value,
                pdb_file=os.path.abspath(tmp_file_path),
                center_x=job.center_x,
                center_y=job.center_y,
                center_z=job.center_z,
                size_x=job.size_x,
                size_y=job.size_y,
                size_z=job.size_z,
                batch_size=job.batch_size,
                minscore=job.minscore,
                epochs=job.epochs
            )

            job.change_sampling_size(5)

            job.save(cascade=True)

            self.jobs.append(job)

    async def jobs_setup_errors(self) -> List[JobValidationError]:
        jobs_errors = []

        for job in self.jobs:
            errors = job.input_errors()
            if errors:
                jobs_errors.append(
                    JobValidationError(
                        job_id=job.id,
                        msg=', '.join([err.message for err in errors])
                    )
                )

        return jobs_errors

    @property
    def _input_parameter_type(self) -> Type[SmallMoleculesDesignLearningInput]:
        return SmallMoleculesDesignLearningInput

    @property
    def _output_parameter_type(self) -> Type[SmallMoleculesDesignLearningOutput]:
        return SmallMoleculesDesignLearningOutput
