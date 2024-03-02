from biobuddy.api_models import SendMessageToBioBuddyResponse, SendMessageToBioBuddyRequest

from biobuddy.loggers import logger

from openai import OpenAI

__all__ = ['send_message']

client = OpenAI()


def send_message(request: SendMessageToBioBuddyRequest) -> SendMessageToBioBuddyResponse:
    completion = client.chat.completions.create(
        model="gpt-3.5-turbo",
        messages=[
            {"role": "system",
             "content": "You are a biobuddy. You are a research assistant that helps creating new drugs, researching biology, pulling information from different sources and running experiments."},
            {"role": "user", "content": request.message_content}
        ],
        tools=request.tools
    )

    return SendMessageToBioBuddyResponse(chatgpt_reply=completion.choices[0].message)
