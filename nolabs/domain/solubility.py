import pydantic.dataclasses as pcdataclasses


@pcdataclasses.dataclass
class SolubilityProbability:
    value: float
