from pydantic.dataclasses import dataclass

from nolabs.api_models.drug_discovery import TargetMetaData

@dataclass
class TargetId:
    value: str

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            raise ValueError()

        return self.value == other.value

    def __str__(self):
        return self.value

    def __hash__(self):
        return self.value.__hash__()

@dataclass
class TargetData:
    fasta_content: str
    pdb_content: str | None = None
    pocket_ids: int | None = None

@dataclass
class Target:
    id: str
    metaData: TargetMetaData
    data: TargetData

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            raise ValueError()

        return self.id == other.id

    def __str__(self):
        return (self.id, self.name)

    def __hash__(self):
        return self.id.__hash__()