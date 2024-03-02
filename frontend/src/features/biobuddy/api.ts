import {
  CancelablePromise, BiobuddyService, LoadConversationResponse, SendMessageResponse,
} from 'src/api/client';

export function loadConversationApi(experimentId: string): CancelablePromise<LoadConversationResponse> {
  return BiobuddyService.loadConversationApiV1BiobuddyLoadConversationGet(experimentId);
}

export function sendMessageApi(experimentId: string, message: string): CancelablePromise<SendMessageResponse> {
  return BiobuddyService.sendMessageApiV1BiobuddySendMessagePost(experimentId, message);
}
