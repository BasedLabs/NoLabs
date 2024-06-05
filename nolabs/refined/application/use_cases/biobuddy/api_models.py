from typing import List, Optional, Any
from uuid import UUID

from pydantic import BaseModel, model_validator
from pydantic.dataclasses import dataclass


@dataclass
class CheckBioBuddyEnabledResponse:
    enabled: bool


