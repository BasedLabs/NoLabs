from typing import Optional

import pydantic


@pydantic.dataclasses.dataclass
class ExperimentProperties:
    contig: str
    number_of_designs: int
    input_file_name: str
    input_file_content: str
    timesteps: Optional[int] = None
    hotspots: Optional[str] = None
