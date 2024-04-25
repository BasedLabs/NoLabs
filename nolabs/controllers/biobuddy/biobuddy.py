from typing import Annotated

from fastapi import APIRouter, Depends

from nolabs.api_models.biobuddy import LoadConversationResponse, LoadConversationRequest, CreateMessageResponse, \
    CreateMessageRequest, CheckBioBuddyEnabledResponse, SendQueryResponse, SendQueryRequest, EditMessageResponse, \
    EditMessageRequest
from nolabs.controllers.biobuddy.dependencies import load_conversation_dependency, \
    create_message_dependency, send_query_dependency, check_biobuddy_enabled_dependency, edit_message_dependency
from nolabs.modules.biobuddy.check_biobuddy_enabled_feature import CheckBioBuddyEnabledFeature
from nolabs.modules.biobuddy.edit_message_feature import EditMessageFeature
from nolabs.modules.biobuddy.load_conversation_feature import LoadConversationFeature
from nolabs.modules.biobuddy.save_message_feature import CreateMessageFeature
from nolabs.modules.biobuddy.send_query_feature import SendQueryFeature

router = APIRouter(
    prefix='/api/v1/biobuddy',
    tags=['biobuddy']
)

@router.get('/check_biobuddy_enabled')
async def check_biobuddy_enabled(
        feature: Annotated[
    CheckBioBuddyEnabledFeature, Depends(check_biobuddy_enabled_dependency)]) -> CheckBioBuddyEnabledResponse:
    return feature.handle()

@router.get('/load-conversation')
async def load_conversation(experiment_id: str,
        feature: Annotated[
    LoadConversationFeature, Depends(load_conversation_dependency)]) -> LoadConversationResponse:
    return feature.handle(LoadConversationRequest(experiment_id))


@router.post('/message/create')
async def save_message(experiment_id: str,
                       message_content: str,
        feature: Annotated[
    CreateMessageFeature, Depends(create_message_dependency)]) -> CreateMessageResponse:
    return feature.handle(CreateMessageRequest(experiment_id, message_content))

@router.post('/message/edit')
async def edit_message(experiment_id: str,
                       message_id: str,
                       message_content: str,
        feature: Annotated[
    EditMessageFeature, Depends(edit_message_dependency)]) -> EditMessageResponse:
    return feature.handle(EditMessageRequest(experiment_id, message_id, message_content))

@router.post('/query')
async def send_query(experiment_id: str,
                       message_content: str,
        feature: Annotated[
    SendQueryFeature, Depends(send_query_dependency)]) -> SendQueryResponse:
    return feature.handle(SendQueryRequest(experiment_id, message_content))

