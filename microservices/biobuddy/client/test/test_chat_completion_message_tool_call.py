# coding: utf-8

"""
    Bio Buddy

    No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)

    The version of the OpenAPI document: 0.1.0
    Generated by OpenAPI Generator (https://openapi-generator.tech)

    Do not edit the class manually.
"""  # noqa: E501


import unittest

from biobuddy_microservice.models.chat_completion_message_tool_call import (
    ChatCompletionMessageToolCall,
)


class TestChatCompletionMessageToolCall(unittest.TestCase):
    """ChatCompletionMessageToolCall unit test stubs"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def make_instance(self, include_optional) -> ChatCompletionMessageToolCall:
        """Test ChatCompletionMessageToolCall
        include_option is a boolean, when False only required
        params are included, when True both required and
        optional params are included"""
        # uncomment below to create an instance of `ChatCompletionMessageToolCall`
        """
        model = ChatCompletionMessageToolCall()
        if include_optional:
            return ChatCompletionMessageToolCall(
                id = '',
                function = { },
                type = None
            )
        else:
            return ChatCompletionMessageToolCall(
                id = '',
                function = { },
                type = None,
        )
        """

    def testChatCompletionMessageToolCall(self):
        """Test ChatCompletionMessageToolCall"""
        # inst_req_only = self.make_instance(include_optional=False)
        # inst_req_and_optional = self.make_instance(include_optional=True)


if __name__ == "__main__":
    unittest.main()
