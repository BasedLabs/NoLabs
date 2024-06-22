__all__ = [
    'CheckBioBuddyEnabledFeature',
    'EditMessageFeature',
    'LoadConversationFeature',
    'CreateMessageFeature',
    'SendQueryFeature',
]

import ast
import json
import os
import uuid
from typing import List, Dict
from uuid import UUID

import biobuddy_microservice

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.application.use_cases.biobuddy.api_models import CheckBioBuddyEnabledResponse, EditMessageRequest, \
    EditMessageResponse, LoadConversationRequest, LoadConversationResponse, CreateMessageRequest, CreateMessageResponse, \
    SendQueryRequest, SendQueryResponse, GetAvailableFunctionCallsResponse
from nolabs.application.use_cases.biobuddy.api_models import Message as ApiMessage
from nolabs.application.use_cases.biobuddy.api_models import RegularMessage as ApiRegularMessage
from nolabs.application.use_cases.biobuddy.api_models import FunctionCall as ApiFunctionCall
from nolabs.application.use_cases.biobuddy.api_models import FunctionParam as ApiFunctionParam
from nolabs.domain.models.biobuddy import Chat, TextMessage, FunctionCall, UserRoleEnum, FunctionParam, \
    FunctionCallMessage


class CheckBioBuddyEnabledFeature:
    def __init__(self):
        pass

    def handle(self) -> CheckBioBuddyEnabledResponse:
        enable_biobuddy = os.getenv('ENABLE_BIOBUDDY', 'false').lower() == 'true'
        return CheckBioBuddyEnabledResponse(enabled=enable_biobuddy)


class EditMessageFeature:
    def __init__(self):
        pass

    def handle(self, request: EditMessageRequest) -> EditMessageResponse:
        assert request

        experiment_id = request.experiment_id
        message_id = request.message_id

        chat = Chat.objects(experiment_id=experiment_id).first()
        if not chat:
            raise NoLabsException(ErrorCodes.invalid_experiment_id, 'Experiment not found')

        chat.edit_message(message_id, request.message_content)
        chat.save()

        return EditMessageResponse(ApiMessage(id=message_id,
                                              role='user',
                                              message=ApiRegularMessage(
                                                  content=request.message_content
                                              ),
                                              type='text')
                                   )


class LoadConversationFeature:
    def __init__(self):
        pass

    def handle(self, request: LoadConversationRequest) -> LoadConversationResponse:
        assert request

        experiment_id = request.experiment_id

        chat = Chat.objects(experiment_id=experiment_id).first()
        if not chat:
            return LoadConversationResponse(messages=[])

        messages = chat.messages
        api_messages = []
        for message in messages:
            if isinstance(message, TextMessage):
                api_message = ApiMessage(
                    id=message.message_id,
                    role=message.sender,
                    message=ApiRegularMessage(content=message.content),
                    type='text'
                )
            elif isinstance(message, FunctionCallMessage):
                function_calls = [
                    ApiFunctionCall(
                        function_name=function_call.function_name,
                        arguments=[ApiFunctionParam(name=param.name, value=param.value) for param in function_call.arguments]
                    )
                    for function_call in message.function_calls
                ]
                api_message = ApiMessage(
                    id=message.message_id,
                    role=message.sender,
                    message=function_calls,
                    type='function'
                )
            else:
                continue
            api_messages.append(api_message)

        return LoadConversationResponse(messages=api_messages)

class CreateMessageFeature:
    def __init__(self):
        pass

    def handle(self, request: CreateMessageRequest) -> CreateMessageResponse:
        assert request

        experiment_id = request.experiment_id

        message_id = uuid.uuid4()

        chat = Chat.objects(experiment_id=experiment_id).first()
        if not chat:
            chat = Chat(experiment_id=experiment_id, messages=[])

        if request.role == 'user':
            chat.add_text_message(message_id, UserRoleEnum.user, request.message_content)
        else:
            chat.add_text_message(message_id, UserRoleEnum.biobuddy, request.message_content)

        chat.save()

        return CreateMessageResponse(ApiMessage(id=message_id,
                                                role=request.role,
                                                message=ApiRegularMessage(
                                                    content=request.message_content
                                                ),
                                                type='text'))

class GetAvailableFunctionCallsFeature:
    def __init__(self, functions: List):
        self._functions = {function.name: function for function in functions}
        self._tools = self.construct_tools_object()

    def handle(self) -> GetAvailableFunctionCallsResponse:
        return GetAvailableFunctionCallsResponse(function_calls=self._tools)


    def construct_tools_object(self):
        tools = []
        if not self._functions:
            return tools
        for func in self._functions.values():
            tool = {
                "type": "function",
                "function": {
                    "name": func.name,
                    "description": func.description,
                    "parameters": {
                        "type": "object",
                        "properties": {},
                        "required": [],
                    },
                }
            }
            properties = tool["function"]["parameters"]["properties"]
            required_params = tool["function"]["parameters"]["required"]

            for param in func.parameters:
                param_info = {
                    "type": param.type,
                    "description": param.description
                }
                if param.type == "array":
                    assert param.items_type, (f"items_type must "
                                              f"be specified for array parameters, missing in {param.name}")
                    param_info["items"] = {"type": param.items_type}

                properties[param.name] = param_info

                if param.required:
                    required_params.append(param.name)

            tools.append(tool)
        return tools

class SendQueryFeature:
    def __init__(self, biobuddy_microservice: biobuddy_microservice.DefaultApi, functions: List):
        self._biobuddy_microservice = biobuddy_microservice
        self._functions = {function.name: function for function in functions}
        self._tools = self.construct_tools_object()

    def handle(self, request: SendQueryRequest) -> SendQueryResponse:
        assert request

        experiment_id = request.experiment_id

        previous_messages = self.get_previous_messages(experiment_id)

        request = biobuddy_microservice.SendMessageToBioBuddyRequest(message_content=request.query,
                                                                     previous_messages=previous_messages,
                                                                     tools=self._tools)

        # try:
        assistant_message = self._biobuddy_microservice.predict_send_message_post(request)

        if assistant_message.reply_type == "function":
            return self._process_ai_function_calls(experiment_id, assistant_message)
        # except Exception as e:
        #    raise NoLabsException(ErrorCodes.biobuddy_error_generating_response)

        return self._process_regular_ai_response(experiment_id, assistant_message)

    def _process_ai_function_calls(self, experiment_id: UUID, assistant_message) -> SendQueryResponse:
        function_calls = ast.literal_eval(assistant_message.content)
        api_function_calls = []
        regular_function_calls = []
        for function_call in function_calls:
            function_dict = ast.literal_eval(function_call)['function_call']
            function_name = function_dict['name']
            arguments = json.loads(function_dict['arguments'])
            if function_name in self._functions:
                function_call = self._functions[function_name].execute(arguments=arguments)

                regular_function_calls.append(FunctionCall(
                                                           function_name=function_call.function_name,
                                                           arguments=[FunctionParam(name=param.name,
                                                                                      value=str(param.value)) for param in
                                                                     function_call.arguments]))

                api_function_calls.append(ApiFunctionCall(function_name=function_call.function_name,
                                                          arguments=[ApiFunctionParam(name=param.name,
                                                                                      value=param.value) for param in
                                                                     function_call.arguments],
                                                          data=function_call.data))

        message_id = uuid.uuid4()
        chat = Chat.objects(experiment_id=experiment_id).first()
        if not chat:
            chat = Chat(experiment_id=experiment_id, messages=[])

        chat.add_function_call_message(message_id, UserRoleEnum.biobuddy, regular_function_calls)
        chat.save()

        return SendQueryResponse(biobuddy_response=ApiMessage(id=message_id,
                                                              role='assistant',
                                                              message=api_function_calls,
                                                              type='function'))

    def _process_regular_ai_response(self, experiment_id: UUID, assistant_message) -> SendQueryResponse:
        message_id = uuid.uuid4()
        chat = Chat.objects(experiment_id=experiment_id).first()
        if not chat:
            chat = Chat(experiment_id=experiment_id, messages=[])

        chat.add_text_message(message_id, UserRoleEnum.biobuddy, str(assistant_message.content))
        chat.save()

        return SendQueryResponse(biobuddy_response=ApiMessage(id=message_id,
                                                              role=UserRoleEnum.biobuddy.value,
                                                              message=ApiRegularMessage(
                                                                  content=str(assistant_message.content)
                                                              ),
                                                              type='text'))

    def construct_tools_object(self):
        tools = []
        if not self._functions:
            return tools
        for func in self._functions.values():
            tool = {
                "type": "function",
                "function": {
                    "name": func.name,
                    "description": func.description,
                    "parameters": {
                        "type": "object",
                        "properties": {},
                        "required": [],
                    },
                }
            }
            properties = tool["function"]["parameters"]["properties"]
            required_params = tool["function"]["parameters"]["required"]

            for param in func.parameters:
                param_info = {
                    "type": param.type,
                    "description": param.description
                }
                if param.type == "array":
                    assert param.items_type, (f"items_type must "
                                              f"be specified for array parameters, missing in {param.name}")
                    param_info["items"] = {"type": param.items_type}

                properties[param.name] = param_info

                if param.required:
                    required_params.append(param.name)

            tools.append(tool)
        return tools

    def get_previous_messages(self, experiment_id: UUID) -> List[Dict[str, str]]:
        chat = Chat.objects(experiment_id=experiment_id).first()
        if not chat:
            return []

        previous_messages = []
        for message in chat.messages:
            if isinstance(message, TextMessage):
                previous_messages.append({'role': message.sender.value,
                                          'content': message.content})
            elif isinstance(message, FunctionCall):
                functions = [{"name": message.function_name, "arguments": message.arguments}]
                previous_messages.append({'role': UserRoleEnum.biobuddy.value,
                                          'content': f"I called {functions}"})

        return previous_messages
