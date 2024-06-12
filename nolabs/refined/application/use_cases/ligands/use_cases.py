__all__ = [
    'SearchLigandsFeature',
    'GetLigandFeature',
    'DeleteLigandFeature',
    'UploadLigandFeature'
]

import base64
import pathlib
from typing import List, Optional
from uuid import UUID

from mongoengine import Q

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.ligands.api_models import LigandResponse, LigandSearchQuery, \
    UploadLigandRequest, UpdateLigandRequest
from nolabs.refined.domain.models.common import Experiment, Ligand, LigandName


def map_ligand_to_response(ligand: Ligand) -> LigandResponse:
    return LigandResponse(
        id=ligand.iid.value,
        experiment_id=ligand.experiment.iid.value,
        name=str(ligand.name),
        smiles_content=ligand.get_smiles(),
        sdf_content=ligand.get_sdf(),
        link=ligand.link,
        image=base64.b64encode(ligand.image).decode("utf-8"),
        drug_likeness=ligand.drug_likeness,
        designed_ligand_score=ligand.designed_ligand_score
    )


class SearchLigandsFeature:
    async def handle(self, query: LigandSearchQuery) -> List[LigandResponse]:
        db_query = Q()

        if query.name:
            db_query = db_query & Q(name__icontains=query.name)

        if query.experiment_id:
            experiment = Experiment.objects.with_id(query.experiment_id)
            db_query = db_query & Q(experiment=experiment)

        return [map_ligand_to_response(l) for l in Ligand.objects(db_query)]


class GetLigandFeature:
    async def handle(self, ligand_id: UUID) -> Optional[LigandResponse]:
        ligand = Ligand.objects.with_id(ligand_id)

        if not ligand:
            return None

        return map_ligand_to_response(ligand)


class UploadLigandFeature:
    async def handle(self, request: UploadLigandRequest) -> LigandResponse:
        experiment = Experiment.objects.with_id(request.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        if not request.smiles and not request.sdf:
            raise NoLabsException(ErrorCodes.invalid_ligand_content, 'You must specify either smiles or sdf')

        name = request.name

        if request.smiles:
            ext = pathlib.Path(request.smiles.filename).suffix
            if ext not in ['.smiles', '.smi']:
                raise NoLabsException(ErrorCodes.smiles_file_is_invalid)
            name = name or request.smiles.filename

        if request.sdf:
            ext = pathlib.Path(request.sdf.filename).suffix
            if ext not in ['.sdf']:
                raise NoLabsException(ErrorCodes.sdf_file_is_invalid)
            name = name or request.sdf.filename

        ligand_name = LigandName(name)

        ligand = Ligand.create(
            name=ligand_name,
            experiment=experiment,
            smiles_content=await request.smiles.read() if request.smiles else None,
            sdf_content=await request.sdf.read() if request.sdf else None
        )
        ligand.save(cascade=True)

        return map_ligand_to_response(ligand)


class UpdateLigandFeature:
    async def handle(self, request: UpdateLigandRequest) -> LigandResponse:

        ligand = Ligand.objects.with_id(request.ligand_id)

        if not ligand:
            raise NoLabsException(ErrorCodes.ligand_not_found)

        if request.smiles:
            ext = pathlib.Path(request.smiles.filename).suffix
            if ext not in ['.smiles', '.smi']:
                raise NoLabsException(ErrorCodes.smiles_file_is_invalid)
            ligand.set_smiles(await request.smiles.read())

        if request.sdf:
            ext = pathlib.Path(request.sdf.filename).suffix
            if ext not in ['.sdf']:
                raise NoLabsException(ErrorCodes.sdf_file_is_invalid)
            ligand.set_smiles(await request.sdf.read())

        if request.name:
            ligand_name = LigandName(request.name or request.smiles.filename or request.sdf.filename)
            ligand.set_name(name=ligand_name)

        ligand.save(cascade=True)

        return map_ligand_to_response(ligand)


class DeleteLigandFeature:
    async def handle(self, ligand_id: UUID):
        ligand = Ligand.objects.with_id(ligand_id)

        if not ligand:
            return

        ligand.delete()



