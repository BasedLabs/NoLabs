import dataclasses
import json
from typing import List, Dict

from nolabs.api_models.biobuddy import SendMessageRequest, SendMessageResponse
from nolabs.api_models.biobuddy import Message as ApiMessage
from nolabs.domain.experiment import ExperimentId
from nolabs.features.biobuddy.data_models.message import Message, RegularMessage, FunctionCall, FunctionParam
from nolabs.api_models.biobuddy import RegularMessage as ApiRegularMessage
from nolabs.api_models.biobuddy import FunctionCall as ApiFunctionCall
from nolabs.features.biobuddy.file_management import FileManagement

import biobuddy_microservice
from biobuddy_microservice import DefaultApi
from biobuddy_microservice import ApiClient
from biobuddy_microservice import Configuration
from nolabs.features.biobuddy.functions.base_function import BiobuddyFunction
from nolabs.infrastructure.settings import Settings

class SendMessageFeature:
    def __init__(self, settings: Settings, file_management: FileManagement, functions: List[BiobuddyFunction] = None):
        self._file_management = file_management
        self._settings = settings
        self._functions = {function.name: function for function in functions}
        self._tools = self.construct_tools_object()

    def handle(self, request: SendMessageRequest) -> SendMessageResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)

        configuration = Configuration(
            host=self._settings.biobuddy_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)

            previous_messages = self.get_previous_messages(experiment_id)

            request = biobuddy_microservice.SendMessageToBioBuddyRequest(message_content=request.message_content,
                                                                         previous_messages=previous_messages,
                                                                         tools=self._tools)
            assistant_message = api_instance.predict_send_message_post(
                send_message_to_bio_buddy_request=request).chatgpt_reply

            print(assistant_message)
            if assistant_message.tool_calls and assistant_message.tool_calls.actual_instance:
                function_name = assistant_message.tool_calls.actual_instance[0]['function']['name']
                arguments = json.loads(assistant_message.tool_calls.actual_instance[0]['function']['arguments'])
                if function_name in self._functions:
                    function_call = self._functions[function_name].execute(experiment_id=experiment_id,
                                                                           arguments=arguments)
                    self._file_management.update_conversation(experiment_id,
                                                              Message(role='user',
                                                                      message=RegularMessage(
                                                                          content=request.message_content
                                                                      ),
                                                                      type='text')
                                                              )
                    self._file_management.update_conversation(experiment_id,
                                                              Message(role='assistant',
                                                                      message=FunctionCall(
                                                                          function_name=function_call.function_name,
                                                                          parameters=[FunctionParam(name=param.name,
                                                                                                    value=param.value) for param in
                                                                                      function_call.parameters]
                                                                                           ),
                                                                      type="function")
                                                              )
                    return SendMessageResponse(biobuddy_response=ApiMessage(role='assistant',
                                                                            message=ApiFunctionCall(**dataclasses.asdict(function_call)),
                                                                            type="function")
                                               )

            else:
                self._file_management.update_conversation(experiment_id,
                                                          Message(role='user',
                                                                  message=RegularMessage(
                                                                      content=request.message_content
                                                                  ),
                                                                  type='text')
                                                          )
                self._file_management.update_conversation(experiment_id,
                                                          Message(role='assistant',
                                                                  message=RegularMessage(
                                                                      content=str(
                                                                          assistant_message.content.actual_instance)
                                                                  ),
                                                                  type='text')
                                                          )
                return SendMessageResponse(biobuddy_response=ApiMessage(role='assistant',
                                                                        message=ApiRegularMessage(
                                                                            content=str(assistant_message.content.actual_instance)
                                                                        ),
                                                                        type='text')
                                           )

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
            if message.type == "text":
                previous_messages.append({'role': message.role,
                                          'content': str(message.message.content)})
            elif message.type == "function":
                previous_messages.append({'role': message.role,
                                          'content': f"I called {message.message.function_name} with parameters:"
                                                     f"{message.message.parameters}"})

        return previous_messages

