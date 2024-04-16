from leaf import FileObject
from pydantic.dataclasses import dataclass

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.seedwork.domain.value_objects import ValueObject, ReadableValueObject, TReadableContent, GenericUUID


@dataclass
class ProteinName(ValueObject):
    value: str

    _min_length = 1
    _max_length = 1000

    def __init__(self, value: str):
        if not value:
            raise NoLabsException('Protein name cannot be empty', ErrorCodes.invalid_protein_name)

        if len(value) < self._min_length or len(value) > self._max_length:
            raise NoLabsException('Length of protein name must be between 1 and 1000', ErrorCodes.invalid_protein_name)

        self.value = value


@dataclass
class AminoAcidName(ValueObject):
    value: str

    _min_length = 1
    _max_length = 1000

    def __init__(self, value: str):
        if not value:
            raise NoLabsException('Amino acid name cannot be empty', ErrorCodes.invalid_aa_name)

        if len(value) < self._min_length or len(value) > self._max_length:
            raise NoLabsException('Length of amino acid name must be between 1 and 1000', ErrorCodes.invalid_aa_name)

        self.value = value


@dataclass
class ExperimentName(ValueObject):
    value: str

    _min_length = 1
    _max_length = 100

    def __init__(self, value: str):
        if not value:
            NoLabsException.throw(ErrorCodes.invalid_experiment_name)

        if len(value) < self._min_length or len(value) > self._max_length:
            NoLabsException.throw(ErrorCodes.invalid_experiment_name)

        self.value = value


@dataclass
class AminoAcidContent(ReadableValueObject[str]):
    _file: FileObject

    def __init__(self, file: FileObject):
        if not file:
            NoLabsException.throw(ErrorCodes.invalid_aa_content)

        self._file = file

    def read(self) -> TReadableContent:
        return self._file.read_string()

    def write(self, content: TReadableContent):
        self._file.write_string(content)


@dataclass
class ProteinContent(ReadableValueObject[str]):
    _file: FileObject

    def __init__(self, file: FileObject):
        if not file:
            NoLabsException.throw(ErrorCodes.invalid_protein_content)

        self._file = file

    def read(self) -> TReadableContent:
        return self._file.read_string()

    def write(self, content: TReadableContent):
        self._file.write_string(content)


@dataclass
class ProteinId(ValueObject, GenericUUID):
    _value: str

    def __init__(self, uuid: str):
        super().__init__()

        if not self.validate(uuid):
            NoLabsException.throw(ErrorCodes.invalid_protein_id)

        self._value = uuid


@dataclass
class ExperimentId(ValueObject, GenericUUID):
    _value: str

    def __init__(self, uuid: str):
        super().__init__()

        if not self.validate(uuid):
            NoLabsException.throw(ErrorCodes.invalid_experiment_id)

        self._value = uuid


@dataclass
class AminoAcidId(ValueObject, GenericUUID):
    _value: str

    def __init__(self, uuid: str):
        super().__init__()

        if not self.validate(uuid):
            NoLabsException.throw(ErrorCodes.invalid_aa_id)

        self._value = uuid
