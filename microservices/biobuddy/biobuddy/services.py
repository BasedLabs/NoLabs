from typing import Dict, AsyncGenerator, List, Any

from langchain.callbacks import AsyncIteratorCallbackHandler
from langchain_core.output_parsers import JsonOutputParser
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder

from biobuddy.api_models import SendMessageToBioBuddyResponse, SendMessageToBioBuddyRequest
from biobuddy.prompts import generate_strategy_prompt, generate_system_prompt, generate_workflow_prompt
from langchain_core.messages import HumanMessage, SystemMessage, AIMessage
from langchain_openai import ChatOpenAI

from biobuddy.gpt_researcher import get_report
import asyncio

chat_model = ChatOpenAI()
chat_model.model_name = "gpt-4-turbo"
chat_model.temperature = 0.1
chat_model.streaming = True


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

    if "<WORKFLOW>" in completion.content:
        new_strategy_prompt = generate_workflow_prompt(request.message_content)
        completion = runnable.invoke({
            "input": new_strategy_prompt,
            "history": history_messages
        })
        response_content = completion.content
        return SendMessageToBioBuddyResponse(
            reply_type="workflow_reply",
            content=response_content
        )

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


async def send_message_async(request: SendMessageToBioBuddyRequest, stop_tokens: Dict[str, bool]) -> AsyncGenerator[SendMessageToBioBuddyResponse, None]:
    system_message_content = generate_system_prompt()
    history_messages = [SystemMessage(content=system_message_content)]
    for msg in request.previous_messages:
        if msg['role'] == 'user':
            history_messages.append(HumanMessage(content=msg['content']))
        elif msg['role'] == 'assistant':
            history_messages.append(AIMessage(content=msg['content']))

    tools_description = " ".join([f"{(tool['function']['name'], tool['function']['description'])}, " for tool in request.tools])
    strategy_prompt = generate_strategy_prompt(tools_description, request.message_content)

    handler = AsyncIteratorCallbackHandler()
    response_buffer = ""

    # Start the streaming
    async def chat_with_model():
        await chat_model.agenerate(
            messages=[[SystemMessage(content=system_message_content), HumanMessage(content=strategy_prompt)]],
            stop=["<STOP>"],
            callbacks=[handler]
        )

    # Launch the chat model invocation
    asyncio.create_task(chat_with_model())

    # Yield responses as they come in
    async for response in handler.aiter():
        if stop_tokens.get(request.experiment_id):
            break
        response_buffer += response

        # Ensure other coroutines get a chance to run
        await asyncio.sleep(0)

        if "<ACTION>" in response_buffer and "<END_ACTION>" in response_buffer:
            start_idx = response_buffer.index("<ACTION>") + len("<ACTION>")
            end_idx = response_buffer.index("<END_ACTION>")
            action_text = response_buffer[start_idx:end_idx].strip()
            response_buffer = ""
            action_response = await invoke_action(action_text, request.tools)
            yield SendMessageToBioBuddyResponse(reply_type="function", content=action_response)

        # Stream intermediate responses if no complete action is found
        if not ("<ACTION>" in response_buffer and "<END_ACTION>" in response_buffer):
            yield SendMessageToBioBuddyResponse(reply_type="stream", content=response)

    # Handle the final response separately
    if not stop_tokens.get(request.experiment_id):
        if response_buffer:
            yield SendMessageToBioBuddyResponse(reply_type="stream", content=response_buffer)
        yield SendMessageToBioBuddyResponse(reply_type="final", content="<STOP>")

async def invoke_action(action_text: str, tools: List[Dict[str, Any]]) -> str:
    # Use ainvoke to handle the action
    action_response = await chat_model.ainvoke(
        input=[HumanMessage(content=action_text)],
        stop=["<STOP>"],
        functions=[tool['function'] for tool in tools]
    )

    print("ACTION TEXT, ACTION RESPONSE: ", (action_text, action_response))

    return str(action_response.additional_kwargs)
