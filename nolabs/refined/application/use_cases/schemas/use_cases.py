from typing import List

from nolabs.refined.application.use_cases.schemas.api_models import ApiSchema


class GetApiSchemasFeature:
    async def handle(self) -> List[ApiSchema]:
        pass