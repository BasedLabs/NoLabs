from pydantic.dataclasses import dataclass


@dataclass
class File:
    name: str
    content: str


@dataclass
class Message:
    role: str
    content: str | File
    type: str
