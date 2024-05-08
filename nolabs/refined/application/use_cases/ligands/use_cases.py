__all__ = [
    'SearchLigandsFeature',
    'GetLigandFeature',
    'DeleteLigandFeature',
    'UploadLigandFeature'
]

from typing import List, Optional
from uuid import UUID

from mongoengine import Q

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.ligands.api_models import LigandResponse, LigandSearchQuery, \
    UploadLigandRequest
from nolabs.refined.domain.models.common import Experiment, Ligand, LigandName


def map_ligand_to_response(ligand: Ligand) -> LigandResponse:
    return LigandResponse(
        id=ligand.iid.value,
        experiment_id=ligand.experiment.iid.value,
        name=str(ligand.name),
        smiles_content=ligand.get_smiles(),
        sdf_content=ligand.get_sdf(),
        drug_likeness=ligand.drug_likeness,
        designed_ligand_score=ligand.designed_ligand_score
    )


class SearchLigandsFeature:
    async def handle(self, query: LigandSearchQuery) -> List[LigandResponse]:
        db_query = Q()

        if query.name:
            db_query = db_query or Q(name=query.name)

        if query.experiment_id:
            experiment = Experiment.objects.with_id(query.experiment_id)
            db_query = db_query or Q(experiment=experiment)

        return [map_ligand_to_response(l) for l in Ligand.objects(db_query)]


class GetLigandFeature:
    async def handle(self, ligand_id: UUID) -> Optional[LigandResponse]:
        ligand = Ligand.objects.with_id(ligand_id)

        if not ligand:
            return None

        return map_ligand_to_response(ligand)


class UploadLigandFeature:
    async def handle(self, request: UploadLigandRequest) -> LigandResponse:
        ligand_name = LigandName(request.name or request.smiles.filename or request.sdf.filename)

        experiment = Experiment.objects.with_id(request.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        ligand = Ligand.create(
            name=ligand_name,
            experiment=experiment,
            smiles_content=await request.smiles.read() if request.smiles else None,
            sdf_content=await request.sdf.read() if request.sdf else None
        )
        ligand.save(cascade=True)

        return map_ligand_to_response(ligand)


class DeleteLigandFeature:
    async def handle(self, ligand_id: UUID):
        ligand = Ligand.objects.with_id(ligand_id)

        if not ligand:
            return

        ligand.delete()



