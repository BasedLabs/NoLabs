from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Depends

from nolabs.application.biobuddy.api_models import (
    CheckBioBuddyEnabledResponse,
    CreateFunctionCallMessageRequest,
    CreateFunctionCallMessageResponse,
    CreateMessageRequest,
    CreateMessageResponse,
    EditMessageRequest,
    EditMessageResponse,
    FunctionCall,
    GetAvailableFunctionCallsResponse,
    LoadConversationRequest,
    LoadConversationResponse,
    SendQueryRequest,
    SendQueryResponse,
)
from nolabs.application.biobuddy.di import BiobuddyDependencies
from nolabs.application.biobuddy.use_cases import (
    CheckBioBuddyEnabledFeature,
    CreateFunctionCallMessageFeature,
    CreateMessageFeature,
    EditMessageFeature,
    GetAvailableFunctionCallsFeature,
    LoadConversationFeature,
    SendActionQueryFeature,
)

router = APIRouter(prefix="/api/v1/biobuddy", tags=["biobuddy"])


@router.get("/check_biobuddy_enabled")
async def check_biobuddy_enabled(
    feature: Annotated[
        CheckBioBuddyEnabledFeature,
        Depends(BiobuddyDependencies.check_biobuddy_enabled),
    ]
) -> CheckBioBuddyEnabledResponse:
    return feature.handle()


@router.get("/load-conversation")
async def load_conversation(
    experiment_id: UUID,
    feature: Annotated[
        LoadConversationFeature, Depends(BiobuddyDependencies.load_conversation)
    ],
) -> LoadConversationResponse:
    return feature.handle(LoadConversationRequest(experiment_id))


@router.post("/message/create")
async def create_message(
    experiment_id: UUID,
    message_id: UUID,
    message_content: str,
    role: str,
    feature: Annotated[
        CreateMessageFeature, Depends(BiobuddyDependencies.create_message)
    ],
) -> CreateMessageResponse:
    return feature.handle(
        CreateMessageRequest(experiment_id, message_id, message_content, role)
    )


@router.post("/function/create")
async def create_function_call_message(
    experiment_id: UUID,
    message_id: UUID,
    function_call: FunctionCall,
    role: str,
    feature: Annotated[
        CreateFunctionCallMessageFeature,
        Depends(BiobuddyDependencies.create_function_call_message),
    ],
) -> CreateFunctionCallMessageResponse:
    return feature.handle(
        CreateFunctionCallMessageRequest(experiment_id, message_id, function_call, role)
    )


@router.post("/message/edit")
async def edit_message(
    experiment_id: UUID,
    message_id: UUID,
    message_content: str,
    feature: Annotated[EditMessageFeature, Depends(BiobuddyDependencies.edit_message)],
) -> EditMessageResponse:
    return feature.handle(
        EditMessageRequest(experiment_id, message_id, message_content)
    )


@router.post("/query")
async def send_query(
    experiment_id: UUID,
    action_text: str,
    feature: Annotated[
        SendActionQueryFeature, Depends(BiobuddyDependencies.send_query)
    ],
) -> SendQueryResponse:
    return feature.handle(SendQueryRequest(experiment_id, action_text))


@router.get("/tools")
async def get_tools(
    feature: Annotated[
        GetAvailableFunctionCallsFeature, Depends(BiobuddyDependencies.get_tools)
    ]
) -> GetAvailableFunctionCallsResponse:
    return feature.handle()
