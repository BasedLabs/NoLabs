from pydantic.dataclasses import dataclass
from nolabs.refined.domain.common.entities import AminoAcidId


@dataclass(frozen=True)
class AminoAcidCreatedEvent:
    amino_acid_id: AminoAcidId
