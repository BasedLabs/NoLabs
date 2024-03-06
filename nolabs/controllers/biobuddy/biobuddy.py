from typing import Annotated

from fastapi import APIRouter, Depends

from nolabs.api_models.biobuddy import LoadConversationResponse, LoadConversationRequest, SendMessageResponse, \
    SendMessageRequest
from nolabs.controllers.biobuddy.dependencies import load_conversation_dependency, send_message_to_drug_discovery_dependency
from nolabs.features.biobuddy.load_conversation_feature import LoadConversationFeature
from nolabs.features.biobuddy.send_message_feature import SendMessageFeature

router = APIRouter(
    prefix='/api/v1/biobuddy',
    tags=['biobuddy']
)

@router.get('/load-conversation')
async def load_conversation(experiment_id: str,
        feature: Annotated[
    LoadConversationFeature, Depends(load_conversation_dependency)]) -> LoadConversationResponse:
    return feature.handle(LoadConversationRequest(experiment_id))


@router.post('/send-message')
async def send_message(experiment_id: str,
                       message_content: str,
        feature: Annotated[
    SendMessageFeature, Depends(send_message_to_drug_discovery_dependency)]) -> SendMessageResponse:
    return feature.handle(SendMessageRequest(experiment_id, message_content))