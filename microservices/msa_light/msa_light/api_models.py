from __future__ import annotations

from typing import Optional

import dataclasses

@dataclasses.dataclass
class RunMsaPredictionRequest:
    api_url: str
    fasta_contents: str
    job_id: str = None

@dataclasses.dataclass
class RunMsaPredictionResponse:
    msa_contents: str = None
