import {
  CancelablePromise,
  BiobuddyService,
  LoadConversationResponse,
  CreateMessageResponse,
  CheckBioBuddyEnabledResponse,
  SendQueryResponse, EditMessageResponse,
  GetAvailableFunctionCallsResponse,
} from 'src/refinedApi/client';

export function checkBioBuddyEnabled(): CancelablePromise<CheckBioBuddyEnabledResponse> {
  return BiobuddyService.checkBiobuddyEnabledApiV1BiobuddyCheckBiobuddyEnabledGet();
}

export function loadConversationApi(experimentId: string): CancelablePromise<LoadConversationResponse> {
  return BiobuddyService.loadConversationApiV1BiobuddyLoadConversationGet(experimentId);
}

export function saveMessageApi(experimentId: string, message: string, role: string): CancelablePromise<CreateMessageResponse> {
  return BiobuddyService.createMessageApiV1BiobuddyMessageCreatePost(experimentId, message, role);
}

export function editMessageApi(experimentId: string, messageId: string, message: string): CancelablePromise<EditMessageResponse> {
  return BiobuddyService.editMessageApiV1BiobuddyMessageEditPost(experimentId, messageId, message);
}


export function sendQueryApi(experimentId: string, query: string): CancelablePromise<SendQueryResponse> {
  return BiobuddyService.sendQueryApiV1BiobuddyQueryPost(experimentId, query);
}

export function getToolsApi(): CancelablePromise<GetAvailableFunctionCallsResponse> {
  return BiobuddyService.getToolsApiV1BiobuddyToolsGet();
}