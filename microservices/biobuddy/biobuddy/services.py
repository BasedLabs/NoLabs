from langchain_core.output_parsers import JsonOutputParser
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder

from biobuddy.api_models import SendMessageToBioBuddyResponse, SendMessageToBioBuddyRequest
from langchain_core.messages import HumanMessage, SystemMessage, AIMessage
from langchain_openai import ChatOpenAI

from biobuddy.gpt_researcher import get_report
import asyncio

chat_model = ChatOpenAI()
chat_model.model_name = "gpt-4-turbo"
chat_model.temperature = 0.1


def send_message(request: SendMessageToBioBuddyRequest) -> SendMessageToBioBuddyResponse:
    system_message_content = (
        "You are a biobuddy. You are a research assistant that helps creating new drugs, researching biology, "
        "pulling information from different sources and running experiments. "
        "If a user asks you to pull the targets only do so for proteins and use rcsb pdb. If a user asks for "
        "drugs/ligands then use chembl."
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

    tools_description = " ".join([f"{(tool['function']['name'], tool['function']['description'])}, " for tool in request.tools])

    openai_functions = [tool['function'] for tool in request.tools]

    strategy_prompt = (f"Given the available tools ({tools_description}), decide whether a direct reply is sufficient "
                       f"or if a specific plan involving these tools is necessary. If plan is not needed, provide the direct reply. "
                       f"If plan is needed, reply with the list (with a tag <PLAN> before the list) of function calls in the format <PLAN>: ['Call function_1 in order to do X', 'Call function_2 to achieve Y']. "
                       f"If you need to consult external database of knowledge to make a well-structured answer, write only tag <RESEARCH>. "
                       f"\n\nQuery: \"{request.message_content}\"\n\n")


    tool_prompt = f"""You are an assistant that has access to the following set of tools. Here are the names and descriptions for each tool:

    {request.tools}

    Given the user input, return the name and input of the tool to use. Return your response as a JSON blob with 'name' and 'arguments' keys."""

    prompt = ChatPromptTemplate.from_messages(
        [("system", tool_prompt),
         ("human", "{input}")]
    )

    tool_chain = prompt | chat_model | JsonOutputParser()

    completion = runnable.invoke({
        "input": strategy_prompt,
        "history": history_messages
    })

    if "<PLAN>" in completion.content:
        plan_text = completion.content.split("<PLAN>:")[1].strip()
        plan_actions = eval(plan_text)

        function_calls = []
        for action in plan_actions:
            print(action)
            action_completion = chat_model.invoke( [HumanMessage(content=action)], functions=openai_functions)
            function_calls.append(str(action_completion.additional_kwargs))

        response_content = str(function_calls)
        return SendMessageToBioBuddyResponse(
            reply_type="function_calls",
            content=response_content
        )
    elif "<RESEARCH>" in completion.content:
        response_content = asyncio.run(get_report(request.message_content))
        return SendMessageToBioBuddyResponse(
            reply_type="regular_reply", # TODO: change to research_reply and make frontend display it differently (it's markdown)
            content=response_content
        )
    else:
        response_content = completion.content
        return SendMessageToBioBuddyResponse(
            reply_type="regular_reply",
            content=response_content
        )


