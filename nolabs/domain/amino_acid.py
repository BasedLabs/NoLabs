from pydantic.dataclasses import dataclass


@dataclass
class AminoAcid:
    name: str
    sequence: str


