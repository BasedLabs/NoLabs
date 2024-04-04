from biobuddy.api_models import SendMessageToBioBuddyResponse, SendMessageToBioBuddyRequest

from langchain_core.messages import HumanMessage, SystemMessage, AIMessage
from langchain.chains import create_history_aware_retriever
from langchain_openai import ChatOpenAI

chat_model = ChatOpenAI()
chat_model.model_name = "gpt-4"
chat_model.streaming = False
chat_model.temperature = 0.1

__all__ = ['send_message']


def send_message(request: SendMessageToBioBuddyRequest) -> SendMessageToBioBuddyResponse:
    # Initialize the messages list with the SystemMessage
    messages = [SystemMessage(
        content="You are a biobuddy. You are a research assistant that helps creating new drugs, researching biology, "
                "pulling information from different sources and running experiments."
                "If a user asks you to pull the targets only do so for proteins and use rcsb pdb. If a user asks for "
                "drugs/ligands them use chembl. If a user asks for latest research about the illness then use pubmed.")]

    # Add previous messages based on their role
    for msg in request.previous_messages:
        print(msg)
        if msg["role"] == "user":
            messages.append(HumanMessage(content=msg["message"]))
        elif msg["role"] == "assistant":
            messages.append(AIMessage(content=msg["message"]))

    completion = chat_model.invoke(

        input=request.message_content
    )

    return SendMessageToBioBuddyResponse(chatgpt_reply=completion.content)
