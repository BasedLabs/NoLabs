from __future__ import annotations

import dataclasses

@dataclasses.dataclass
class RunMsaPredictionRequest:
    api_url: str
    fasta_contents: str

@dataclasses.dataclass
class RunMsaPredictionResponse:
    msa_contents: str
