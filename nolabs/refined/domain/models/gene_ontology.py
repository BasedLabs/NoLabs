__all__ = [
    'GeneOntologyJob'
]

from typing import List, Dict, Any
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, DictField, EmbeddedDocument, UUIDField, \
    EmbeddedDocumentListField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein


class GeneOntologyJobResult(EmbeddedDocument):
    protein_id: UUID = UUIDField(required=True)
    gene_ontology: Dict[str, Any] = DictField(required=False)


class GeneOntologyJob(Job):
    proteins: List[Protein] = ListField(ReferenceField(Protein, required=False, reverse_delete_rule=PULL))
    gene_ontologies: List[GeneOntologyJobResult] = EmbeddedDocumentListField(GeneOntologyJobResult)

    def set_proteins(self, proteins: List[Protein]):
        if not proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        self.proteins = proteins

    def clear_result(self):
        self.gene_ontologies = []

    def set_result(self, protein: Protein, gene_ontology: Dict[str, Any]):
        if not gene_ontology:
            raise NoLabsException(ErrorCodes.invalid_gene_ontology)

        if not self.proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if protein not in self.proteins:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        if not protein.fasta_content:
            raise NoLabsException(ErrorCodes.protein_amino_acid_sequence_not_found)

        existing_result = [res for res in self.gene_ontologies if res.protein_id == protein.id]
        if existing_result:
            gene_ontology_job_result = existing_result[0]
            gene_ontology_job_result.gene_ontology = gene_ontology
        else:
            gene_ontology_job_result = GeneOntologyJobResult(
                protein_id=protein.id,
                gene_ontology=gene_ontology
            )
            self.gene_ontologies.append(gene_ontology_job_result)
