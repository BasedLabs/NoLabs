from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.refined.application.use_cases.biobuddy.api_models import CheckBioBuddyEnabledResponse, \
    LoadConversationRequest, LoadConversationResponse, CreateMessageResponse, CreateMessageRequest, EditMessageResponse, \
    EditMessageRequest, SendQueryResponse, SendQueryRequest
from nolabs.refined.application.use_cases.biobuddy.di import BiobuddyDependencies
from nolabs.refined.application.use_cases.biobuddy.use_cases import CheckBioBuddyEnabledFeature, \
    LoadConversationFeature, CreateMessageFeature, EditMessageFeature, SendQueryFeature

router = APIRouter(
    prefix='/api/v1/biobuddy',
    tags=['biobuddy']
)

@router.get('/check_biobuddy_enabled')
async def check_biobuddy_enabled(
        feature: Annotated[
    CheckBioBuddyEnabledFeature, Depends(BiobuddyDependencies.check_biobuddy_enabled)]) -> CheckBioBuddyEnabledResponse:
    return feature.handle()

@router.get('/load-conversation')
async def load_conversation(experiment_id: UUID,
        feature: Annotated[
    LoadConversationFeature, Depends(BiobuddyDependencies.load_conversation)]) -> LoadConversationResponse:
    return feature.handle(LoadConversationRequest(experiment_id))


@router.post('/message/create')
async def create_message(experiment_id: UUID,
                       message_content: str,
        feature: Annotated[
    CreateMessageFeature, Depends(BiobuddyDependencies.create_message)]) -> CreateMessageResponse:
    return feature.handle(CreateMessageRequest(experiment_id, message_content))

@router.post('/message/edit')
async def edit_message(experiment_id: UUID,
                       message_id: UUID,
                       message_content: str,
        feature: Annotated[
    EditMessageFeature, Depends(BiobuddyDependencies.edit_message)]) -> EditMessageResponse:
    return feature.handle(EditMessageRequest(experiment_id, message_id, message_content))

@router.post('/query')
async def send_query(experiment_id: UUID,
                       message_content: str,
        feature: Annotated[
    SendQueryFeature, Depends(BiobuddyDependencies.send_query)]) -> SendQueryResponse:
    return feature.handle(SendQueryRequest(experiment_id, message_content))