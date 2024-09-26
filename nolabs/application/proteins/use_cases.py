__all__ = [
    "SearchProteinsContentFeature",
    "GetProteinFeature",
    "DeleteProteinFeature",
    "UploadProteinFeature",
    "UpdateProteinFeature",
]

import pathlib
from typing import List, Optional
from uuid import UUID

from mongoengine import Q

from nolabs.application.proteins.api_models import (
    ProteinContentResponse,
    ProteinLocalisationResponse,
    ProteinMetadataResponse,
    ProteinSearchMetadataQuery,
    ProteinSearchQuery,
    UpdateProteinRequest,
    UploadProteinRequest,
    UploadProteinResponse,
)
from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import Experiment, Protein, ProteinName


def map_protein_to_response(protein: Protein) -> ProteinContentResponse:
    return ProteinContentResponse(
        id=protein.iid.value,
        experiment_id=protein.experiment.iid.value,
        name=str(protein.name),
        fasta_content=protein.get_fasta(),
        pdb_content=protein.get_pdb(),
        binding_pockets=protein.binding_pockets,
        localisation=(
            ProteinLocalisationResponse(
                cytosolic=protein.localisation.cytosolic,
                mitochondrial=protein.localisation.mitochondrial,
                nuclear=protein.localisation.nuclear,
                other=protein.localisation.other,
                extracellular=protein.localisation.extracellular,
            )
            if protein.localisation
            else None
        ),
        gene_ontology=protein.gene_ontology,
        soluble_probability=protein.soluble_probability,
        msa=protein.get_msa(),
        md_pdb_content=protein.get_md(),
        fasta_name=protein.name.fasta_name,
        pdb_name=protein.name.pdb_name,
        link=str(protein.link),
    )


def map_to_protein_metadata_response(protein: Protein) -> ProteinMetadataResponse:
    return ProteinMetadataResponse(
        id=protein.iid.value,
        experiment_id=protein.experiment.iid.value,
        name=str(protein.name),
        binding_pockets=protein.binding_pockets,
        localisation=(
            ProteinLocalisationResponse(
                cytosolic=protein.localisation.cytosolic,
                mitochondrial=protein.localisation.mitochondrial,
                nuclear=protein.localisation.nuclear,
                other=protein.localisation.other,
                extracellular=protein.localisation.extracellular,
            )
            if protein.localisation
            else None
        ),
        gene_ontology=protein.gene_ontology,
        soluble_probability=protein.soluble_probability,
        fasta_name=protein.name.fasta_name,
        pdb_name=protein.name.pdb_name,
        link=str(protein.link),
    )


def map_to_upload_protein_response(protein: Protein) -> UploadProteinResponse:
    return UploadProteinResponse(
        id=protein.iid.value,
        experiment_id=protein.experiment.iid.value,
        name=str(protein.name),
        binding_pockets=protein.binding_pockets,
        localisation=(
            ProteinLocalisationResponse(
                cytosolic=protein.localisation.cytosolic,
                mitochondrial=protein.localisation.mitochondrial,
                nuclear=protein.localisation.nuclear,
                other=protein.localisation.other,
                extracellular=protein.localisation.extracellular,
            )
            if protein.localisation
            else None
        ),
        gene_ontology=protein.gene_ontology,
        soluble_probability=protein.soluble_probability,
        fasta_name=protein.name.fasta_name,
        pdb_name=protein.name.pdb_name,
        link=str(protein.link),
    )


class SearchProteinsContentFeature:
    async def handle(self, query: ProteinSearchQuery) -> List[ProteinContentResponse]:
        db_query = Q()

        if (
            not query.all
            and not query.name
            and not query.experiment_id
            and not query.ids
        ):
            return []

        if query.all:
            return [map_protein_to_response(p) for p in Protein.objects()]

        if query.ids:
            proteins = Protein.objects(id__in=query.ids)

            return [map_protein_to_response(p) for p in proteins]

        if query.name:
            db_query = db_query & Q(name__icontains=query.name)

        if query.experiment_id:
            experiment = Experiment.objects.with_id(query.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            db_query = db_query & Q(experiment=experiment)

        return [map_protein_to_response(p) for p in Protein.objects(db_query)]


class SearchProteinsMetadataFeature:
    async def handle(
        self, query: ProteinSearchMetadataQuery
    ) -> List[ProteinMetadataResponse]:
        db_query = Q()

        if (
            not query.all
            and not query.name
            and not query.experiment_id
            and not query.ids
        ):
            return []

        if query.all:
            return [map_to_protein_metadata_response(p) for p in Protein.objects()]

        if query.ids:
            proteins = Protein.objects(id__in=query.ids)

            return [map_to_protein_metadata_response(p) for p in proteins]

        if query.name:
            db_query = db_query & Q(name__icontains=query.name)

        if query.experiment_id:
            experiment = Experiment.objects.with_id(query.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            db_query = db_query & Q(experiment=experiment)

        return [map_to_protein_metadata_response(p) for p in Protein.objects(db_query)]


class GetProteinFeature:
    async def handle(self, protein_id: UUID) -> Optional[ProteinContentResponse]:
        protein = Protein.objects.with_id(protein_id)

        if not protein:
            return None

        return map_protein_to_response(protein)


class GetProteinMetadataFeature:
    async def handle(self, protein_id: UUID) -> Optional[ProteinMetadataResponse]:
        protein = Protein.objects.with_id(protein_id)

        if not protein:
            return None

        return map_to_protein_metadata_response(protein)


class UploadProteinFeature:
    async def handle(self, request: UploadProteinRequest) -> UploadProteinResponse:
        experiment = Experiment.objects.with_id(request.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        if not request.fasta and not request.pdb:
            raise NoLabsException(
                ErrorCodes.invalid_protein_content,
                "You must specify either pdb or fasta",
            )

        if request.fasta and not request.fasta.filename.endswith(".fasta"):
            raise NoLabsException(ErrorCodes.pdb_file_is_invalid)

        if request.pdb and not request.pdb.filename.endswith(".pdb"):
            raise NoLabsException(ErrorCodes.fasta_file_is_invalid)

        name = request.name

        if request.fasta:
            name = request.fasta.filename or name

        if request.pdb:
            name = request.pdb.filename or name

        protein_name = ProteinName(name)

        protein = Protein.create(
            name=protein_name,
            experiment=experiment,
            fasta_content=await request.fasta.read() if request.fasta else None,
            pdb_content=await request.pdb.read() if request.pdb else None,
        )
        protein.save(cascade=True)

        return map_to_upload_protein_response(protein)


class UpdateProteinFeature:
    async def handle(self, request: UpdateProteinRequest) -> ProteinContentResponse:

        protein = Protein.objects.with_id(request.protein_id)

        if not protein:
            raise NoLabsException(ErrorCodes.protein_not_found)

        if request.fasta:
            if pathlib.Path(request.fasta.filename).suffix != ".fasta":
                raise NoLabsException(ErrorCodes.fasta_file_is_invalid)
            protein.set_fasta(await request.fasta.read())

        if request.pdb:
            if pathlib.Path(request.pdb.filename).suffix != ".pdb":
                raise NoLabsException(ErrorCodes.pdb_file_is_invalid)
            protein.set_pdb(await request.pdb.read())

        if request.name:
            protein_name = ProteinName(
                request.name or request.pdb.filename or request.fasta.filename
            )
            protein.set_name(name=protein_name)

        protein.save(cascade=True)

        return map_protein_to_response(protein)


class DeleteProteinFeature:
    async def handle(self, protein_id: UUID):
        protein = Protein.objects.with_id(protein_id)

        if not protein:
            return

        protein.delete()
