from __future__ import annotations

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
import uuid
from functools import lru_cache
from pathlib import Path
from typing import Union, Dict, Any, List
from uuid import UUID
from Bio import SeqIO
import io

from mongoengine import DateTimeField, Document, ReferenceField, CASCADE, EmbeddedDocument, \
    FloatField, EmbeddedDocumentField, BinaryField, UUIDField, DictField, ListField, IntField, \
    Q, StringField
from pydantic import model_validator
from pydantic.dataclasses import dataclass
from rdkit import Chem
from typing_extensions import Self
from nolabs.utils.generate_2d_drug import generate_png_from_smiles

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
        try:
            uuid.UUID(str(self.value))
        except ValueError:
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

        if len(self.value) > 1000:
            raise NoLabsException(ErrorCodes.invalid_experiment_name)
        return self


class Experiment(Document, Entity):
    id: UUID = UUIDField(primary_key=True)
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

    def __hash__(self):
        return self.iid.__hash__()

    def __eq__(self, other):
        if isinstance(other, Experiment):
            return self.iid == other.iid

        return False


@dataclass
class ProteinName(ValueObjectString):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_protein_name)

        if len(self.value) > 1000:
            raise NoLabsException(ErrorCodes.invalid_protein_name)

        value = Path(self.value).stem

        self.value = value

        return self

    @property
    def fasta_name(self):
        return self.value + '.fasta'

    @property
    def pdb_name(self):
        return self.value + '.pdb'


@dataclass
class ProteinId(ValueObjectUUID):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        try:
            uuid.UUID(str(self.value))
        except ValueError:
            raise NoLabsException(ErrorCodes.invalid_protein_id)

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
    id: UUID = UUIDField(primary_key=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name: ProteinName = ValueObjectStringField(required=True, factory=ProteinName)
    fasta_content: bytes | None = BinaryField(required=False)
    pdb_content: bytes | None = BinaryField(required=False)
    localisation: LocalisationProbability | None = EmbeddedDocumentField(LocalisationProbability, required=False)
    gene_ontology: Dict[str, Any] | None = DictField(required=False)
    soluble_probability: SolubleProbability | None = ValueObjectFloatField(required=False, factory=SolubleProbability)
    msa: bytes | None = BinaryField(required=False)

    binding_pockets: List[int] = ListField(IntField())
    md_content: bytes | None = BinaryField(required=False)

    source_binding_protein = ReferenceField('Protein', required=False)
    binding_ligand: Ligand = ReferenceField('Ligand', required=False)
    sdf_content: bytes | None = BinaryField(required=False)
    minimized_affinity: float | None = FloatField(required=False)
    scored_affinity: float | None = FloatField(required=False)
    confidence: float | None = FloatField(required=False)
    plddt_array: List[int] = ListField(IntField, required=False)

    '''
    Conformations content
    '''

    def get_msa(self) -> str | None:
        if self.msa:
            return self.msa.decode('utf-8')

        return None

    def get_protein_binders(self) -> List['ProteinBinder']:
        return ProteinBinder.objects(Q(protein1=self) or Q(protein2=self))

    def get_ligand_binders(self) -> List['LigandBinder']:
        return LigandBinder.objects(protein=self)

    @property
    def iid(self) -> ProteinId:
        return ProteinId(self.id)

    def set_name(self, name: ProteinName):
        if not name:
            raise NoLabsException(ErrorCodes.invalid_protein_name)

        self.name = name

    def set_fasta(self, fasta_content: Union[bytes, str]):
        if not fasta_content:
            raise NoLabsException(ErrorCodes.protein_fasta_is_empty)

        if isinstance(fasta_content, str):
            fasta_content = fasta_content.encode('utf-8')

        self.fasta_content = fasta_content

    def get_fasta(self) -> str | None:
        if self.fasta_content:
            return self.fasta_content.decode('utf-8')

        return None

    def get_amino_acid_sequence(self) -> str | None:
        fasta = self.get_fasta()
        if fasta:
            res = ''
            for chain in SeqIO.parse(io.StringIO(fasta), 'fasta'):
                res += str(chain.seq)

            return res

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

    @classmethod
    def create(cls, experiment: Experiment,
               name: ProteinName,
               fasta_content: Union[bytes, str, None] = None,
               pdb_content: Union[bytes, str, None] = None,
               *args,
               **kwargs):
        if not id:
            raise NoLabsException(ErrorCodes.invalid_protein_id)

        if not name:
            raise NoLabsException(ErrorCodes.invalid_protein_name)
        if not experiment:
            raise NoLabsException(ErrorCodes.invalid_experiment_id)

        if fasta_content and isinstance(fasta_content, str):
            fasta_content = fasta_content.encode('utf-8')

        if pdb_content and isinstance(fasta_content, str):
            pdb_content = pdb_content.encode('utf-8')

        if not fasta_content and not pdb_content:
            raise NoLabsException(ErrorCodes.protein_initialization_error)

        proteins = cls.objects(name=name.value, experiment=experiment)
        if proteins:
            protein: Protein = proteins[0]

            if fasta_content and not protein.fasta_content:
                protein.set_fasta(fasta_content)

            if pdb_content and not protein.pdb_content:
                protein.set_pdb(pdb_content)

            protein.set_name(name)
            return protein

        protein = Protein(
            id=ProteinId(uuid.uuid4()).value,
            experiment=experiment,
            name=name,
            *args,
            **kwargs
        )

        if fasta_content:
            protein.set_fasta(fasta_content)

        if pdb_content:
            protein.set_pdb(pdb_content)

        EventDispatcher.raise_event(ProteinCreatedEvent(protein))

        return protein

    def set_localisation_probability(self, localisation: LocalisationProbability):
        if not localisation:
            raise NoLabsException(ErrorCodes.invalid_localisation_probability)

        self.localisation = localisation

    def add_protein_binder(self, protein: 'Protein'):
        if self == protein:
            raise NoLabsException(ErrorCodes.protein_cannot_be_binder_to_itself)

        for binder in self.get_protein_binders():
            if binder.protein1 == protein or binder.protein2 == protein:
                return

        binder = ProteinBinder(protein1=self, protein2=protein)
        binder.save()

    def set_msa(self, msa: bytes | bool):
        if not msa:
            raise NoLabsException(ErrorCodes.invalid_msa)

        if isinstance(msa, str):
            msa = msa.encode('utf-8')

        self.msa = msa

    def __hash__(self):
        return self.iid.__hash__

    def __eq__(self, other):
        if not isinstance(other, Protein):
            return False

        return self.iid == other.iid


@dataclass
class LigandId(ValueObjectUUID):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        try:
            uuid.UUID(str(self.value))
        except ValueError:
            raise NoLabsException(ErrorCodes.invalid_ligand_id)

        return self


@dataclass
class LigandName(ValueObjectString):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_ligand_name)

        value = Path(self.value).stem

        self.value = value

        return self


@dataclass
class DrugLikenessScore(ValueObjectFloat):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_drug_likeness_score)

        if self.value < 0 or self.value > 1.0:
            raise NoLabsException(ErrorCodes.invalid_drug_likeness_score,
                                  'Drug likeness score must be in a range [0,1.0]')

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
            raise NoLabsException(ErrorCodes.invalid_designed_ligand_score,
                                  'Designed ligand score must be in a range [0,1.0]')

        return self


class Ligand(Document, Entity):
    id = UUIDField(db_field='_id', primary_key=True, required=True)
    experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name = ValueObjectStringField(required=False, factory=LigandName)
    smiles_content = BinaryField(required=False)
    sdf_content = BinaryField(required=False)
    drug_likeness = FloatField(required=False)
    designed_ligand_score = FloatField(required=False)
    link = StringField(required=False)  # New field for link
    image = BinaryField(required=False)  # New field for image

    def __hash__(self):
        return self.iid.__hash__()

    def __eq__(self, other):
        if isinstance(other, Ligand):
            return self.iid == other.iid
        return False

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
        self.image = self._generate_image(smiles_content)  # Generate image

    def _generate_image(self, smiles_content: bytes) -> bytes:
        image = generate_png_from_smiles(smiles_content.decode('utf-8'))
        from io import BytesIO
        buffer = BytesIO()
        image.save(buffer, format="PNG")
        return buffer.getvalue()


    def get_sdf(self) -> str | None:
        if self.sdf_content:
            return self.sdf_content.decode('utf-8')
        return None

    def get_smiles(self) -> str | None:
        if self.smiles_content:
            return self.smiles_content.decode('utf-8')
        return None

    def set_sdf(self, sdf: Union[str, bytes]):
        if isinstance(sdf, str):
            sdf = sdf.encode('utf-8')
        self.sdf_content = sdf
        self.set_smiles_from_sdf(sdf)  # Set smiles from sdf

    def set_smiles_from_sdf(self, sdf: Union[str, bytes]):
        if isinstance(sdf, bytes):
            sdf = sdf.decode('utf-8')
        mol = Chem.MolFromMolBlock(sdf)
        if mol is None:
            raise NoLabsException(ErrorCodes.sdf_file_is_invalid, 'Invalid SDF content')
        smiles = Chem.MolToSmiles(mol)
        self.set_smiles(smiles)

    def set_drug_likeness_score(self, score: DrugLikenessScore):
        if not score:
            raise NoLabsException(ErrorCodes.invalid_drug_likeness_score)
        self.drug_likeness = score

    def set_designed_ligand_score(self, score: DesignedLigandScore):
        if not score:
            raise NoLabsException(ErrorCodes.invalid_designed_ligand_score)
        self.designed_ligand_score = score

    @classmethod
    def create(cls, experiment: Experiment,
               name: LigandName | None = None,
               smiles_content: Union[bytes, str, None] = None,
               sdf_content: Union[bytes, str, None] = None,
               link: str | None = None,
               *args,
               **kwargs) -> 'Ligand':  # Added link parameter
        if not name:
            raise NoLabsException(ErrorCodes.invalid_ligand_name)
        if not experiment:
            raise NoLabsException(ErrorCodes.invalid_experiment_id)

        if not smiles_content and not sdf_content:
            raise NoLabsException(ErrorCodes.ligand_initialization_error,
                                  'Cannot create a ligand without smiles and sdf content')

        if smiles_content and isinstance(smiles_content, str):
            smiles_content = smiles_content.encode()

        if sdf_content and isinstance(sdf_content, str):
            sdf_content = sdf_content.encode()

        ligands = cls.objects(name=name.value, experiment=experiment)
        if ligands:
            ligand: Ligand = ligands[0]

            if smiles_content:
                ligand.set_smiles(smiles_content)

            if sdf_content:
                ligand.set_sdf(sdf_content)

            ligand.set_name(name)
            ligand.link = link
            return ligand

        if 'id' not in kwargs:
            id = LigandId(uuid.uuid4()).value
        else:
            id = kwargs.get('id')
            if isinstance(id, LigandId):
                id = id.value

        ligand = Ligand(
            id=id,
            experiment=experiment,
            name=name,
            smiles_content=smiles_content,
            sdf_content=sdf_content,
            link=link,  # Set link
            *args,
            **kwargs
        )

        if smiles_content:
            ligand.image = generate_png_from_smiles(smiles_content.decode('utf-8'))  # Generate image
        elif sdf_content:
            ligand.set_smiles_from_sdf(sdf_content.decode('utf-8'))  # Set smiles and image from sdf content

        return ligand

    def add_binding(self,
                    protein: 'Protein',
                    sdf_content: bytes | None = None,
                    minimized_affinity: float | None = None,
                    scored_affinity: float | None = None,
                    confidence: float | None = None,
                    plddt_array: List[int] | None = None,
                    pdb_content: Union[bytes, str, None] = None) -> 'Protein':
        if not plddt_array:
            plddt_array = []

        if not protein:
            raise NoLabsException(ErrorCodes.protein_is_undefined)

        if isinstance(sdf_content, str):
            sdf_content = sdf_content.encode()

        if isinstance(pdb_content, str):
            pdb_content = pdb_content.encode()

        complexes = Protein.objects(ligand=self, protein=protein)

        if complexes:
            complex = complexes[0]
        else:
            complex = Protein.create(
                experiment=self.experiment,
                name=ProteinName(f'{str(protein.name)}-{str(self.name)}-complex'),
                pdb_content=pdb_content,
                fasta_content=protein.fasta_content
            )

        complex.sdf_content = sdf_content
        complex.scored_affinity = scored_affinity if scored_affinity else complex.scored_affinity
        complex.minimized_affinity = minimized_affinity if minimized_affinity else complex.minimized_affinity
        complex.confidence = confidence if confidence else complex.confidence
        complex.plddt_array = plddt_array if plddt_array else complex.plddt_array
        complex.pdb_content = pdb_content if pdb_content else complex.pdb_content

        complex.save()

        return complex

    def get_bindings(self) -> List['LigandBinder']:
        return LigandBinder.objects(ligand=self)


@dataclass
class JobId(ValueObjectUUID):
    value: UUID

    @model_validator(mode='after')
    def post_root(self) -> Self:
        try:
            uuid.UUID(str(self.value))
        except ValueError:
            raise NoLabsException(ErrorCodes.invalid_job_id)

        return self


@dataclass
class JobName(ValueObjectString):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_job_name)

        if len(self.value) > 100:
            raise NoLabsException(ErrorCodes.invalid_job_name)
        return self


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


class ProteinBinder(Document):
    protein1: Protein = ReferenceField(Protein, required=True)
    protein2: Protein = ReferenceField(Protein, required=True)


class LigandBinder(Document):
    protein: Protein = ReferenceField(Protein, required=True, reverse_delete_rule=CASCADE)
    ligand: Ligand = ReferenceField(Ligand, required=True, reverse_delete_rule=CASCADE)
    sdf_content: bytes | None = BinaryField(required=False)
    minimized_affinity: float | None = FloatField(required=False)
    scored_affinity: float | None = FloatField(required=False)
    confidence: float | None = FloatField(required=False)
    plddt_array: List[int] = ListField(IntField, required=False)
    pdb_content: bytes | None = BinaryField(required=False)


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
