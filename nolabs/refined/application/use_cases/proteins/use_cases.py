__all__ = [
    'SearchProteinsFeature',
    'GetProteinFeature',
    'DeleteProteinFeature',
    'UploadProteinFeature'
]

from typing import List, Optional
from uuid import UUID

from mongoengine import Q

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.use_cases.proteins.api_models import ProteinSearchQuery, ProteinResponse, \
    ProteinLocalisationResponse, UploadProteinRequest
from nolabs.refined.domain.models.common import Protein, Experiment, ProteinName


def map_protein_to_response(protein: Protein) -> ProteinResponse:
    return ProteinResponse(
        id=protein.iid.value,
        experiment_id=protein.experiment.iid.value,
        name=str(protein.name),
        fasta_content=protein.get_fasta(),
        pdb_content=protein.get_pdb(),
        localisation=ProteinLocalisationResponse(
            cytosolic=protein.localisation.cytosolic,
            mitochondrial=protein.localisation.mitochondrial,
            nuclear=protein.localisation.nuclear,
            other=protein.localisation.other,
            extracellular=protein.localisation.extracellular
        ) if protein.localisation else None,
        gene_ontology=protein.gene_ontology,
        soluble_probability=protein.soluble_probability,
        msa=protein.get_msa(),
        md_pdb_content=protein.get_md()
    )


class SearchProteinsFeature:
    async def handle(self, query: ProteinSearchQuery) -> List[ProteinResponse]:
        db_query = Q()

        if query.name:
            db_query = db_query or Q(name=query.name)

        if query.experiment_id:
            experiment = Experiment.objects.with_id(query.experiment_id)
            db_query = db_query or Q(experiment=experiment)

        return [map_protein_to_response(p) for p in Protein.objects(db_query)]


class GetProteinFeature:
    async def handle(self, protein_id: UUID) -> Optional[ProteinResponse]:
        protein = Protein.objects.with_id(protein_id)

        if not protein:
            return None

        return map_protein_to_response(protein)


class UploadProteinFeature:
    async def handle(self, request: UploadProteinRequest) -> ProteinResponse:
        protein_name = ProteinName(request.name or request.pdb.filename or request.fasta.filename)

        experiment = Experiment.objects.with_id(request.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        if not request.fasta and not request.pdb:
            raise NoLabsException(ErrorCodes.invalid_protein_content, 'You must specify either pdb or fasta')

        protein = Protein.create(
            name=protein_name,
            experiment=experiment,
            fasta_content=await request.fasta.read() if request.fasta else None,
            pdb_content=await request.pdb.read() if request.pdb else None
        )
        protein.save(cascade=True)

        return map_protein_to_response(protein)


class DeleteProteinFeature:
    async def handle(self, protein_id: UUID):
        protein = Protein.objects.with_id(protein_id)

        if not protein:
            return

        protein.delete()



