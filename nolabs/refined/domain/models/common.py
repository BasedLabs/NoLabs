__all__ = ['Protein',
           'ProteinId',
           'ProteinName',
           'Experiment',
           'Job',
           'JobId',
           'JobName',
           'LocalisationProbability',
           'ProteinCreatedEvent']

import datetime
import os
import uuid
from functools import cached_property, lru_cache
from pathlib import Path
from typing import Union, Dict, Any, List
from uuid import UUID

from mongoengine import DateTimeField, Document, ReferenceField, CASCADE, EmbeddedDocument, \
    FloatField, EmbeddedDocumentField, BinaryField, UUIDField, DictField, ListField, PULL, IntField
from pydantic import model_validator
from pydantic.dataclasses import dataclass
from typing_extensions import Self

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.event_dispatcher import EventDispatcher
from nolabs.refined.infrastructure.mongo_fields import ValueObjectStringField, ValueObjectFloatField
from nolabs.seedwork.domain.entities import Entity
from nolabs.seedwork.domain.events import DomainEvent
from nolabs.seedwork.domain.value_objects import ValueObject, ValueObjectString, ValueObjectUUID, ValueObjectFloat


@dataclass
class ExperimentId(ValueObjectUUID):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_experiment_id)
        return self

    def __eq__(self, other):
        if not isinstance(other, ExperimentId):
            return False

        return self.value == other.value


@dataclass
class ExperimentName(ValueObjectString):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_experiment_name)

        if len(self.value) < 1 or len(self.value) > 1000:
            raise NoLabsException(ErrorCodes.invalid_experiment_name)
        return self


class Experiment(Document, Entity):
    id: UUID = UUIDField(primary_key=True, required=True)
    name: ExperimentName = ValueObjectStringField(required=True, factory=ExperimentName)
    created_at: datetime.datetime = DateTimeField(default=datetime.datetime.utcnow)

    def __init__(self, id: ExperimentId, name: ExperimentName, created_at: datetime.datetime | None = None, *args,
                 **kwargs):
        if not id:
            raise NoLabsException(ErrorCodes.invalid_experiment_id)

        if not name:
            raise NoLabsException(ErrorCodes.invalid_experiment_name)

        created_at = created_at if created_at else datetime.datetime.now(tz=datetime.timezone.utc)

        super().__init__(id=id.value if isinstance(id, ExperimentId) else id, name=name, created_at=created_at, *args,
                         **kwargs)

    @property
    def iid(self) -> ExperimentId:
        return ExperimentId(self.id)

    def set_name(self, name: ExperimentName):
        if not name:
            raise NoLabsException(ErrorCodes.invalid_experiment_name)

        self.name = name


@dataclass
class ProteinName(ValueObjectString):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_protein_name)

        if len(self.value) < 1 or len(self.value) > 1000:
            raise NoLabsException(ErrorCodes.invalid_protein_name)

        value = Path(self.value).stem

        self.value = value

        return self

    @property
    @lru_cache
    def fasta_name(self):
        return self.value + '.fasta'

    @property
    @lru_cache
    def pdb_name(self):
        return self.value + '.pdb'


@dataclass
class ProteinId(ValueObjectUUID):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_protein_id)

        return self


@dataclass
class LigandId(ValueObjectUUID):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_ligand_id)

        return self


@dataclass
class LigandName(ValueObjectString):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_ligand_name)

        return self


@dataclass
class DrugLikenessScore(ValueObjectFloat):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_drug_likeness_score)

        if self.value < 0 or self.value > 1.0:
            raise NoLabsException(ErrorCodes.invalid_drug_likeness_score, 'Drug likeness score must be in a range [0,1.0]')

        return self


@dataclass
class DesignedLigandScore(ValueObjectFloat):
    """
    Average weighted score of a designed ligand.
    """
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_designed_ligand_score)

        if self.value < 0 or self.value > 1.0:
            raise NoLabsException(ErrorCodes.invalid_designed_ligand_score, 'Designed ligand score must be in a range [0,1.0]')

        return self


class Ligand(Document, Entity):
    id: UUID = UUIDField(db_field='_id', primary_key=True, required=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name: LigandName | None = ValueObjectStringField(required=False, factory=LigandName)
    smiles_content: bytes | None = BinaryField(required=True)

    drug_likeness: DrugLikenessScore | None = ValueObjectFloatField(factory=DrugLikenessScore, required=False)
    drug_score: DesignedLigandScore | None = FloatField(factory=DesignedLigandScore, required=False)

    def __init__(
            self,
            id: LigandId,
            experiment: Experiment,
            name: LigandName,
            smiles_content: Union[bytes, str, None] = None,
            *args,
            **kwargs
    ):
        if not id:
            raise NoLabsException(ErrorCodes.invalid_ligand_id)
        if not name:
            raise NoLabsException(ErrorCodes.invalid_ligand_name)
        if not experiment:
            raise NoLabsException(ErrorCodes.invalid_experiment_id)

        if isinstance(smiles_content, str):
            smiles_content = smiles_content.encode('utf-8')

        super().__init__(id=id.value if isinstance(id, LigandId) else id,
                         experiment=experiment,
                         name=name,
                         smiles_content=smiles_content,
                         *args, **kwargs)

        EventDispatcher.raise_event(LigandCreatedEvent(self))

    @property
    def iid(self) -> LigandId:
        return LigandId(self.id)

    def set_name(self, name: LigandName | None):
        self.name = name

    def set_smiles(self, smiles_content: Union[bytes, str]):
        if not smiles_content:
            raise NoLabsException(ErrorCodes.invalid_smiles)

        if isinstance(smiles_content, str):
            smiles_content = smiles_content.encode('utf-8')

        self.smiles_content = smiles_content

    def get_smiles(self) -> str | None:
        if self.smiles_content:
            return self.smiles_content.decode('utf-8')

        return None

    def set_drug_likeness_score(self, score: DrugLikenessScore):
        if not score:
            raise NoLabsException(ErrorCodes.invalid_drug_likeness_score)

        self.drug_likeness = score

    def set_designed_ligand_score(self, score: DesignedLigandScore):
        if not score:
            raise NoLabsException(ErrorCodes.invalid_designed_ligand_score)

        self.drug_score = score

    @classmethod
    def create(cls, experiment: Experiment,
               name: LigandName | None = None,
               smiles_content: Union[bytes, str, None] = None):
        ligands = cls.objects.filter(name=name.value, experiment=experiment)
        if ligands:
            ligand: Ligand = ligands[0]

            if smiles_content:
                ligand.set_smiles(smiles_content)

            ligand.set_name(name)
            return ligand

        return Ligand(
            LigandId(uuid.uuid4()),
            experiment=experiment,
            name=name,
            smiles_content=smiles_content
        )


@dataclass
class JobId(ValueObjectUUID):
    value: UUID

    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_job_id)
        return self


@dataclass
class JobName(ValueObjectString):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_job_name)

        if len(self.value) < 1 or len(self.value) > 100:
            raise NoLabsException(ErrorCodes.invalid_job_name)
        return self


class SolubleProbability(ValueObjectFloat):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_solubility_probability)

        if self.value < 0 or self.value > 1.0:
            raise NoLabsException(ErrorCodes.invalid_solubility_probability)

        return self


class LocalisationProbability(EmbeddedDocument, ValueObject):
    cytosolic: float = FloatField(required=True)
    mitochondrial: float = FloatField(required=True)
    nuclear: float = FloatField(required=True)
    other: float = FloatField(required=True)
    extracellular: float = FloatField(required=True)

    def __init__(self, cytosolic: float, mitochondrial: float,
                 nuclear: float, other: float, extracellular: float, *args, **kwargs):
        values = [
            cytosolic,
            mitochondrial,
            nuclear,
            other,
            extracellular
        ]

        for value in values:
            if not value:
                raise NoLabsException(ErrorCodes.invalid_localisation_probability)
            if value < 0 or value > 1.0:
                raise NoLabsException(ErrorCodes.invalid_localisation_probability)

        super().__init__(
            cytosolic=cytosolic,
            mitochondrial=mitochondrial,
            nuclear=nuclear,
            other=other,
            extracellular=extracellular,
            *args,
            **kwargs)


class Protein(Document, Entity):
    id: UUID = UUIDField(db_field='_id', primary_key=True, required=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name: ProteinName = ValueObjectStringField(required=True, factory=ProteinName)
    fasta_content: bytes | None = BinaryField(required=False)
    pdb_content: bytes | None = BinaryField(required=False)
    localisation: LocalisationProbability | None = EmbeddedDocumentField(LocalisationProbability, required=False)
    gene_ontology: Dict[str, Any] | None = DictField(required=False)
    soluble_probability: SolubleProbability | None = ValueObjectFloatField(required=False, factory=SolubleProbability)
    binding_pockets: List[int] = ListField(IntField(), required=False)
    manual_binding_pockets: List[int] = ListField(IntField(), required=False)
    msa: bytes = BinaryField(required=False)

    binders = ListField(ReferenceField('Protein', required=False, reverse_delete_rule=PULL))
    '''
    Conformations content
    '''
    md_content: bytes | None = BinaryField(required=False)

    def __init__(
            self,
            id: ProteinId,
            experiment: Experiment,
            name: ProteinName,
            fasta_content: Union[bytes, str, None] = None,
            pdb_content: Union[bytes, str, None] = None,
            *args,
            **kwargs
    ):
        if not id:
            raise NoLabsException(ErrorCodes.invalid_protein_id)
        if not name:
            raise NoLabsException(ErrorCodes.invalid_protein_name)
        if not experiment:
            raise NoLabsException(ErrorCodes.invalid_experiment_id)

        if isinstance(fasta_content, str):
            fasta_content = fasta_content.encode('utf-8')

        if isinstance(fasta_content, str):
            pdb_content = pdb_content.encode('utf-8')

        super().__init__(id=id.value if isinstance(id, ProteinId) else id,
                         experiment=experiment,
                         name=name,
                         fasta_content=fasta_content,
                         pdb_content=pdb_content,
                         *args, **kwargs)

        EventDispatcher.raise_event(ProteinCreatedEvent(self))

    @property
    def iid(self) -> ProteinId:
        return ProteinId(self.id)

    def set_name(self, name: ProteinName):
        if not name:
            raise NoLabsException(ErrorCodes.invalid_protein_name)
        self.name = name

    def set_fasta(self, fasta_content: Union[bytes, str]):
        if not fasta_content:
            raise NoLabsException(ErrorCodes.invalid_protein_content)

        if isinstance(fasta_content, str):
            fasta_content = fasta_content.encode('utf-8')

        self.fasta_content = fasta_content

    def get_fasta(self) -> str | None:
        if self.fasta_content:
            return self.fasta_content.decode('utf-8')

        return None

    def set_md(self, md_content: Union[bytes, str]):
        if not md_content:
            raise NoLabsException(ErrorCodes.invalid_protein_content)

        if isinstance(md_content, str):
            md_content = md_content.encode('utf-8')

        self.md_content = md_content

    def get_md(self) -> str | None:
        if self.md_content:
            return self.md_content.decode('utf-8')

        return None

    def set_pdb(self, pdb_content: Union[bytes, str]):
        if not pdb_content:
            raise NoLabsException(ErrorCodes.invalid_protein_content)

        if isinstance(pdb_content, str):
            pdb_content = pdb_content.encode('utf-8')

        self.pdb_content = pdb_content

    def get_pdb(self) -> str | None:
        if self.pdb_content:
            return self.pdb_content.decode('utf-8')

        return None

    def set_gene_ontology(self, gene_ontology: Dict[str, Any]):
        if not gene_ontology:
            raise NoLabsException(ErrorCodes.invalid_gene_ontology)

        self.gene_ontology = gene_ontology

    def set_solubility_probability(self, soluble_probability: SolubleProbability):
        if not soluble_probability:
            raise NoLabsException(ErrorCodes.invalid_solubility_probability)

        self.soluble_probability = soluble_probability

    def set_binding_pockets(self, binding_pockets: List[int]):
        if not binding_pockets:
            raise NoLabsException(ErrorCodes.empty_binding_pockets)

        self.binding_pockets = binding_pockets

    def set_manual_binding_pockets(self, binding_pockets: List[int]):
        if not binding_pockets:
            raise NoLabsException(ErrorCodes.empty_binding_pockets)

        self.manual_binding_pockets = binding_pockets

    @classmethod
    def create(cls, experiment: Experiment,
               name: ProteinName,
               fasta_content: Union[bytes, str, None] = None,
               pdb_content: Union[bytes, str, None] = None):
        proteins = cls.objects.filter(name=name.value, experiment=experiment)
        if proteins:
            protein: Protein = proteins[0]

            if fasta_content:
                protein.set_fasta(fasta_content)

            if pdb_content:
                protein.set_pdb(pdb_content)

            protein.set_name(name)
            return protein

        return Protein(
            ProteinId(uuid.uuid4()),
            experiment=experiment,
            name=name,
            fasta_content=fasta_content,
            pdb_content=pdb_content
        )

    def set_localisation_probability(self, localisation: LocalisationProbability):
        if not localisation:
            raise NoLabsException(ErrorCodes.invalid_localisation_probability)

        self.localisation = localisation

    def set_protein_binder(self, protein: 'Protein'):
        if self == protein:
            raise NoLabsException(ErrorCodes.protein_cannot_be_binder_to_itself)

        for binder in self.binders:
            if binder == protein:
                return

        self.binders.append(protein)

    def set_msa(self, msa: bytes | bool):
        if not msa:
            raise NoLabsException(ErrorCodes.invalid_msa)

        if isinstance(msa, str):
            msa = msa.encode('utf-8')

        self.msa = msa


    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        if not isinstance(other, Protein):
            return False

        return self.id == other.id


class Job(Document, Entity):
    id: UUID = UUIDField(db_field='_id', primary_key=True, required=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name: JobName = ValueObjectStringField(required=True, factory=JobName)
    created_at: datetime.datetime = DateTimeField(default=datetime.datetime.utcnow)

    meta = {
        'allow_inheritance': True
    }

    def __init__(self, id: JobId, name: JobName, experiment: Experiment, *args, **kwargs):
        if not id:
            raise NoLabsException(ErrorCodes.invalid_job_id)
        if not name:
            raise NoLabsException(ErrorCodes.invalid_job_name)
        if not experiment:
            raise NoLabsException(ErrorCodes.invalid_experiment_id)

        super().__init__(id=id.value if isinstance(id, JobId) else id, name=name, experiment=experiment, *args,
                         **kwargs)

    def set_name(self, name: JobName):
        if not name:
            raise NoLabsException(ErrorCodes.invalid_job_name)

        self.name = name

    @property
    def iid(self) -> JobId:
        return JobId(self.id)


# region events
class ProteinCreatedEvent(DomainEvent):
    protein: Protein

    def __init__(self, protein: Protein):
        self.protein = protein


class LigandCreatedEvent(DomainEvent):
    ligand: Ligand

    def __init__(self, ligand: Ligand):
        self.ligand = ligand

# endregion
