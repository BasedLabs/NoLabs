import ast
import dataclasses
import json
from typing import List, Dict

from nolabs.api_models.biobuddy import SendQueryRequest, SendQueryResponse
from nolabs.api_models.biobuddy import Message as ApiMessage
from nolabs.domain.experiment import ExperimentId
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.modules.biobuddy.data_models.message import Message, RegularMessage, FunctionCall, FunctionParam, \
    MessageType
from nolabs.api_models.biobuddy import RegularMessage as ApiRegularMessage
from nolabs.api_models.biobuddy import FunctionCall as ApiFunctionCall
from nolabs.modules.biobuddy.file_management import FileManagement

import biobuddy_microservice
from biobuddy_microservice import DefaultApi
from biobuddy_microservice import ApiClient
from biobuddy_microservice import Configuration
from nolabs.modules.biobuddy.functions.base_function import BiobuddyFunction
from nolabs.infrastructure.settings import Settings
from nolabs.utils import generate_uuid


class SendQueryFeature:
    def __init__(self, settings: Settings, file_management: FileManagement, functions: List[BiobuddyFunction]):
        self._file_management = file_management
        self._settings = settings
        self._functions = {function.name: function for function in functions}
        self._tools = self.construct_tools_object()

    def handle(self, request: SendQueryRequest) -> SendQueryResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)

        configuration = Configuration(
            host=self._settings.biobuddy_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)

            previous_messages = self.get_previous_messages(experiment_id)

            try:
                biobuddy_request = biobuddy_microservice.SendMessageToBioBuddyRequest(
                    message_content=request.query,
                    previous_messages=previous_messages,
                    tools=self._tools)

                assistant_message = api_instance.predict_send_message_post(
                    send_message_to_bio_buddy_request=biobuddy_request)

                print(assistant_message)

                if assistant_message.reply_type == "function":
                    return self._process_ai_function_calls(experiment_id, assistant_message)
            except Exception as e:
                raise NoLabsException(ErrorCodes.biobuddy_error_generating_response)

            return self._process_regular_ai_response(experiment_id, assistant_message)

    def _process_ai_function_calls(self,
                                   experiment_id: ExperimentId,
                                   assistant_message: biobuddy_microservice.SendMessageToBioBuddyResponse
                                   ) -> SendQueryResponse:
        function_calls = ast.literal_eval(assistant_message.content)
        api_function_calls = []
        regular_function_calls = []
        for function_call in function_calls:
            function_dict = ast.literal_eval(function_call)['function_call']
            function_name = function_dict['name']
            arguments = json.loads(function_dict['arguments'])
            if function_name in self._functions:
                function_call = self._functions[function_name].execute(experiment_id=experiment_id,
                                                                       arguments=arguments)

                regular_function_calls.append(FunctionCall(function_name=function_call.function_name,
                                                           parameters=[FunctionParam(name=param.name,
                                                                                     value=param.value) for param in
                                                                       function_call.parameters]
                                                           ))
                api_function_calls.append(ApiFunctionCall(**dataclasses.asdict(function_call)))

        message_id = generate_uuid()
        self._file_management.append_message_to_conversation(experiment_id,
                                                             Message(id=message_id,
                                                                     role='assistant',
                                                                     message=regular_function_calls,
                                                                     type=MessageType.FUNCTIONS)
                                                             )

        return SendQueryResponse(biobuddy_response=ApiMessage(id=message_id,
                                                              role='assistant',
                                                              message=api_function_calls,
                                                              type=MessageType.FUNCTIONS.value))

    def _process_regular_ai_response(self, experiment_id: ExperimentId,
                                     assistant_message: biobuddy_microservice.SendMessageToBioBuddyResponse
                                     ) -> SendQueryResponse:
        message_id = generate_uuid()
        self._file_management.append_message_to_conversation(experiment_id,
                                                             Message(id=message_id,
                                                                     role='assistant',
                                                                     message=RegularMessage(
                                                                         content=str(
                                                                             assistant_message.content)
                                                                     ),
                                                                     type=MessageType.TEXT)
                                                             )
        return SendQueryResponse(biobuddy_response=ApiMessage(id=message_id,
                                                              role='assistant',
                                                              message=ApiRegularMessage(
                                                                  content=str(assistant_message.content)
                                                              ),
                                                              type=MessageType.TEXT.value))

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
                    assert param.items_type, f"items_type must be specified for array parameters, missing in {param.name}"
                    param_info["items"] = {"type": param.items_type}

                properties[param.name] = param_info

                if param.required:
                    required_params.append(param.name)

            tools.append(tool)
        return tools

    def get_previous_messages(self, experiment_id: ExperimentId) -> List[Dict[str, str]]:
        messages = self._file_management.load_conversation(experiment_id)
        previous_messages = []
        for message in messages:
            if message.type == MessageType.TEXT:
                previous_messages.append({'role': message.role,
                                          'content': str(message.message.content)})
            elif message.type == MessageType.FUNCTIONS:
                functions = [(str(idx) + function.function_name, function.parameters) for idx, function in
                             enumerate(message.message)]
                previous_messages.append({'role': message.role,
                                          'content': f"I called {functions}"})

        return previous_messages
