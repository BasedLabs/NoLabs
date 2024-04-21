import uuid
from typing import Optional, Union, Dict, Any, Mapping, List

from bson import Binary
from pymongo.database import Database
from pymongo.collection import Collection

from nolabs.modules.drug_discovery.data_models.result import JobId
from nolabs.refined.domain.common.entities import JobName, JobType, JobState, AminoAcidId, ProteinId, Protein, \
    ProteinName, ProteinContent, AminoAcid, AminoAcidName, AminoAcidContent
from nolabs.refined.domain.common.localisation_entities import LocalisationJob, LocalisationProbability


def _map_to_localisation_job_output(output: Dict[str, Any] | None = None) -> List[LocalisationProbability]:
    return [
        LocalisationProbability(amino_acid_id=AminoAcidId(item['amino_acid_id']),
                                cytosolic_proteins=item['cytosolic_proteins'],
                                mitochondrial_proteins=item['mitochondrial_proteins'],
                                nuclear_proteins=item['nuclear_proteins'],
                                other_proteins=item['other_proteins'],
                                extracellular_secreted_proteins=item['extracellular_secreted_proteins'],
                                )
        for item in output['probabilities']
    ]


def _map_to_localisation_job(job: Mapping[str, Any]) -> LocalisationJob:
    return LocalisationJob(
        id=JobId(job['_id']),
        name=JobName(job['name']),
        job_type=JobType(job['job_type']),
        created_at=job['created_at'],
        state=JobState(job['job_state']),
        amino_acid_ids=[
            AminoAcidId(aa_id) for aa_id in job['amino_acids']],
        output=_map_to_localisation_job_output(job['output'])
    )


class Repository:
    db: Database
    experiments: Collection
    jobs: Collection
    proteins: Collection
    amino_acids: Collection
    files: Collection
    arbitrary: Collection

    def __init__(self, db: Database):
        self.db = db
        self.experiments = db.experiments
        self.jobs = db.jobs
        self.files = db.files
        self.proteins = db.proteins
        self.arbitrary = db.arbitrary

    def get_protein(self, id: ProteinId) -> Protein | None:
        protein_mongo = self.proteins.find_one({'_id': id.value})
        if not protein_mongo:
            return None

        protein_name = ProteinName(protein_mongo['name'])
        content = ProteinContent(protein_mongo['content'])

        protein = Protein(id=id, name=protein_name, content=content)
        protein.clear_events()
        return protein

    def insert_protein(self, protein: Protein):
        self.proteins.insert_one({
            '_id': protein.id.value,
            'name': protein.name.value,
            'content': Binary(protein.content.value)
        })

    def count_amino_acids(self, name: AminoAcidName) -> int:
        return self.amino_acids.count_documents({'name': name.value})

    def get_amino_acid(self, id: AminoAcidId) -> AminoAcid | None:
        aa_mongo = self.amino_acids.find_one({'_id': id.value})
        if not aa_mongo:
            return None

        aa_name = AminoAcidName(aa_mongo['name'])
        content = AminoAcidContent(aa_mongo['content'])

        aa = AminoAcid(id=id, name=aa_name, content=content)
        aa.clear_events()
        return aa

    def get_amino_acids_by_ids(self, ids: List[AminoAcidId]) -> Dict[AminoAcidId, AminoAcid]:
        aa_mongo = self.amino_acids.find({'_id': {'$in': ids}})
        return {
            aa.id: AminoAcid(id=aa.id, name=aa.name, content=aa.content)
            for aa in aa_mongo
        }

    def insert_amino_acid(self, amino_acid: AminoAcid):
        self.amino_acids.insert_one({
            '_id': amino_acid.id.value,
            'name': amino_acid.name.value,
            'content': Binary(amino_acid.content.value)
        })

    def get_localisation_job(self, job_id: JobId) -> LocalisationJob | None:
        job = self.jobs.find_one({'_id': job_id})
        if not job:
            return None
        job = _map_to_localisation_job(job)
        job.clear_events()
        return job

    def insert_arbitrary(self, id: uuid.UUID, d: Dict[str, Any]):
        d['_id'] = id
        self.arbitrary.insert_one(d)

    def get_arbitrary(self, id: uuid.UUID) -> Dict[str, Any] | None:
        return self.arbitrary.find_one({'_id': id})
