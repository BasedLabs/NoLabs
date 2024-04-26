from langchain_core.output_parsers import JsonOutputParser
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder

from biobuddy.api_models import SendMessageToBioBuddyResponse, SendMessageToBioBuddyRequest
from biobuddy.prompts import generate_strategy_prompt, generate_system_prompt
from langchain_core.messages import HumanMessage, SystemMessage, AIMessage
from langchain_openai import ChatOpenAI

from biobuddy.gpt_researcher import get_report
import asyncio

chat_model = ChatOpenAI()
chat_model.model_name = "gpt-4-turbo"
chat_model.temperature = 0.1


def send_message(request: SendMessageToBioBuddyRequest) -> SendMessageToBioBuddyResponse:
    system_message_content = generate_system_prompt()
    history_messages = [SystemMessage(content=system_message_content)]
    for msg in request.previous_messages:
        if msg['role'] == 'user':
            history_messages.append(HumanMessage(content=msg['content']))
        elif msg['role'] == 'assistant':
            history_messages.append(AIMessage(content=msg['content']))

    prompt_template = ChatPromptTemplate.from_messages([
        MessagesPlaceholder(variable_name="history"),
        ("human", "{input}"),
    ])

    runnable = prompt_template | chat_model

    tools_description = " ".join([f"{(tool['function']['name'], tool['function']['description'])}, " for tool in request.tools])

    openai_functions = [tool['function'] for tool in request.tools]

    strategy_prompt = generate_strategy_prompt(tools_description, request.message_content)

    completion = runnable.invoke({
        "input": strategy_prompt,
        "history": history_messages
    })

    if "<PLAN>" in completion.content:
        plan_text = completion.content.split("<PLAN>:")[1].strip()
        plan_actions = eval(plan_text)

        function_calls = []
        for action in plan_actions:
            action_completion = chat_model.invoke( [HumanMessage(content=action)], functions=openai_functions)
            function_calls.append(str(action_completion.additional_kwargs))

        response_content = str(function_calls)
        return SendMessageToBioBuddyResponse(
            reply_type="function",
            content=response_content
        )
    elif "<RESEARCH>" in completion.content:
        response_content = asyncio.run(get_report(request.message_content))
        return SendMessageToBioBuddyResponse(
            reply_type="regular_reply",
            content=response_content
        )
    response_content = completion.content
    return SendMessageToBioBuddyResponse(
        reply_type="regular_reply",
        content=response_content
    )


