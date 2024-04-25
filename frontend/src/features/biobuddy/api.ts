import {
  CancelablePromise,
  BiobuddyService,
  LoadConversationResponse,
  SaveMessageResponse,
  CheckBioBuddyEnabledResponse,
  SendQueryResponse, EditMessageResponse,
} from 'src/api/client';

export function checkBioBuddyEnabled(): CancelablePromise<CheckBioBuddyEnabledResponse> {
  return BiobuddyService.checkBiobuddyEnabledApiV1BiobuddyCheckBiobuddyEnabledGet();
}

export function loadConversationApi(experimentId: string): CancelablePromise<LoadConversationResponse> {
  return BiobuddyService.loadConversationApiV1BiobuddyLoadConversationGet(experimentId);
}

export function saveMessageApi(experimentId: string, message: string): CancelablePromise<SaveMessageResponse> {
  return BiobuddyService.saveMessageApiV1BiobuddyMessageCreatePost(experimentId, message);
}

export function editMessageApi(experimentId: string, messageId: string, message: string): CancelablePromise<EditMessageResponse> {
  return BiobuddyService.editMessageApiV1BiobuddyMessageEditPost(experimentId, messageId, message);
}


export function sendQueryApi(experimentId: string, query: string): CancelablePromise<SendQueryResponse> {
  return BiobuddyService.sendQueryApiV1BiobuddyQueryPost(experimentId, query);
}
