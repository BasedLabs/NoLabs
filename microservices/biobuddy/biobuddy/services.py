from biobuddy.api_models import SendMessageToBioBuddyResponse, SendMessageToBioBuddyRequest
from langchain_core.messages import HumanMessage, SystemMessage, AIMessage
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder

chat_model = ChatOpenAI()
chat_model.model_name = "gpt-4"
chat_model.streaming = False
chat_model.temperature = 0.1

__all__ = ['send_message']


def send_message(request: SendMessageToBioBuddyRequest) -> SendMessageToBioBuddyResponse:
    system_message_content = (
        "You are a biobuddy. You are a research assistant that helps creating new drugs, researching biology, "
        "pulling information from different sources and running experiments. "
        "If a user asks you to pull the targets only do so for proteins and use rcsb pdb. If a user asks for "
        "drugs/ligands then use chembl. If a user asks for latest research about the illness then use pubmed."
    )

    history_messages = [SystemMessage(content=system_message_content)]
    for msg in request.previous_messages:
        if 'user' in msg.keys():
            history_messages.append(HumanMessage(content=msg['user']))
        elif 'assistant' in msg.keys():
            history_messages.append(AIMessage(content=msg['assistant']))

    prompt_template = ChatPromptTemplate.from_messages([
        MessagesPlaceholder(variable_name="history"),
        ("human", "{input}"),
    ])

    runnable = prompt_template | chat_model

    completion = runnable.invoke({
        "input": request.message_content,
        "history": history_messages
    })

    return SendMessageToBioBuddyResponse(chatgpt_reply=completion.content)
