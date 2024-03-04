from biobuddy.api_models import SendMessageToBioBuddyResponse, SendMessageToBioBuddyRequest

from biobuddy.loggers import logger

from openai import OpenAI

__all__ = ['send_message']

client = OpenAI()


def send_message(request: SendMessageToBioBuddyRequest) -> SendMessageToBioBuddyResponse:
    messages = request.previous_messages
    messages.append({"role": "system",
             "content": "You are a biobuddy. You are a research assistant that helps creating new drugs, researching biology, pulling information from different sources and running experiments."
                        "If a user asks you to pull the targets only do so for proteins and use rcsb pdb. If a user asks for drugs/ligands them use chembl. If a user asks for latest research about the illness then use pubmed."})
    messages.append( {"role": "user", "content": request.message_content})
    completion = client.chat.completions.create(
        model="gpt-4",
        messages=messages,
        tools=request.tools
    )

    return SendMessageToBioBuddyResponse(chatgpt_reply=completion.choices[0].message)
