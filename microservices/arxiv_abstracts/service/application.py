import os

from dotenv import load_dotenv

load_dotenv(".env")

from functools import lru_cache
from typing import List, AsyncGenerator
from settings import settings

from api_models import ChatResponse, DocumentContext, MessageType, ChatRequest
from langchain_chroma import Chroma
from langchain_core.documents import Document
from langchain_core.messages import SystemMessage
from langchain_core.tools import tool
from langchain_huggingface import HuggingFaceEmbeddings
from langchain_openai import ChatOpenAI
from langgraph.checkpoint.sqlite.aio import AsyncSqliteSaver
from langgraph.constants import END
from langgraph.prebuilt import ToolNode, tools_condition
from langgraph.graph import StateGraph, MessagesState
from logging import getLogger

logger = getLogger('application')

config = {
    "configurable": {"thread_id": "arxiv"}
}


class State(MessagesState):
    context: List[Document]


@lru_cache
def get_vector_store():
    embeddings = HuggingFaceEmbeddings(model_name="sentence-transformers/all-mpnet-base-v2")
    return Chroma(
        collection_name="arxiv_abstracts",
        persist_directory=settings.chroma_db_path,  # This should be the same directory as before
        embedding_function=embeddings
    )

logger.info('Initialing vector store')

get_vector_store()


@lru_cache
def get_llm():
    return ChatOpenAI(openai_api_key=settings.openai_api_key, temperature=0.2)


@tool(response_format="content_and_artifact")
async def retrieve(query: str):
    """Retrieve information related to a query."""
    vector_store = get_vector_store()
    retrieved_docs = await vector_store.asimilarity_search(query, k=3)
    serialized = "\n\n".join(
        (f"Source: {doc.metadata}\n" f"Content: {doc.page_content}")
        for doc in retrieved_docs
    )
    return serialized, retrieved_docs


async def query_or_respond(state: State):
    llm = get_llm()
    llm_with_tools = llm.bind_tools([retrieve])
    response = await llm_with_tools.ainvoke(state["messages"])
    return {"messages": [response]}


async def generate(state: MessagesState):
    recent_tool_messages = []
    for message in reversed(state["messages"]):
        if message.type == "tool":
            recent_tool_messages.append(message)
        else:
            break
    tool_messages = recent_tool_messages[::-1]

    # Format into prompt
    docs_content = "\n\n".join(doc.content for doc in tool_messages)
    system_message_content = (
        "You are a scientific assistant for question-answering tasks based on documents context. "
        "Use the following arXiv.org abstracts to generate concise summary and answer the question. "
        "If you don't know the answer, say that you don't know and add that you suggest user to rephrase his question so it can better searchable in ChromaDb."
        "Use six sentences maximum and keep the answer concise."
        "\n\n"
        f"{docs_content}"
    )
    conversation_messages = [
        message
        for message in state["messages"]
        if message.type in ("human", "system")
           or (message.type == "ai" and not message.tool_calls)
    ]
    prompt = [SystemMessage(system_message_content)] + conversation_messages

    # Run
    llm = get_llm()
    print(str(prompt))
    response = await llm.ainvoke(prompt)
    context = []
    for tool_message in tool_messages:
        context.extend(tool_message.artifact)
    return {"messages": [response], "context": context}


@lru_cache
def get_graph_builder():
    tools = ToolNode([retrieve])
    graph_builder = StateGraph(MessagesState)

    graph_builder.add_node(query_or_respond)
    graph_builder.add_node(tools)
    graph_builder.add_node(generate)

    graph_builder.set_entry_point("query_or_respond")
    graph_builder.add_conditional_edges(
        "query_or_respond",
        tools_condition,
        {END: END, "tools": "tools"},
    )
    graph_builder.add_edge("tools", "generate")
    graph_builder.add_edge("generate", END)

    return graph_builder


def tool_to_document(tool_data):
    if isinstance(tool_data, dict):
        return DocumentContext(
            page_content=tool_data['page_content'],
            title=tool_data['metadata']['title'],
            id=tool_data['metadata']['id']
        )
    return DocumentContext(
        page_content=tool_data.page_content,
        title=tool_data.metadata['title'],
        id=tool_data.metadata['id']
    )


def to_chat_response(state):
    message_type = MessageType.human
    if state.type == 'human':
        message_type = MessageType.human
    if state.type == 'ai':
        message_type = MessageType.ai
    if state.type == 'tool':
        message_type = MessageType.tool

    if 'artifact' in state.model_fields_set:
        documents = [
            tool_to_document(x) for x in state.artifact
        ]
        return ChatResponse(id=state.id,
                            content=None,
                            context=documents,
                            message_type=message_type)
    content = state.content
    if content:
        return ChatResponse(
            id=state.id,
            content=state.content,
            context=[],
            message_type=message_type
        )
    return None


async def inference(request: ChatRequest) -> AsyncGenerator[ChatResponse, None]:
    graph_builder = get_graph_builder()

    human_question_processed = False

    async with AsyncSqliteSaver.from_conn_string("checkpoint.db") as memory:
        graph = graph_builder.compile(checkpointer=memory)

        async for step in graph.astream(
                {"messages": [{"role": "user", "content": request.message}]},
                stream_mode="values",
                config=config
        ):
            if not human_question_processed:
                human_question_processed = True
                continue
            response = step['messages'][-1]
            response = to_chat_response(response)
            if response:
                yield response


def set_openai_api_key(api_key: str):
    settings.openai_api_key = api_key


async def history() -> List[ChatResponse]:
    graph_builder = get_graph_builder()

    res = []

    async with AsyncSqliteSaver.from_conn_string("checkpoint.db") as memory:
        graph = graph_builder.compile(checkpointer=memory)
        async for holder in graph.aget_state_history(config=config):
            res = []
            for state in holder[0]['messages']:
                response = to_chat_response(state)
                if response:
                    res.append(response)
            return res

    return res

def clear_history():
    os.remove('checkpoint.db')

async def main_history():
    return history()


async def main_inference():
    async for value in inference(request=ChatRequest(message='What is my name?')):
        print(value)


if __name__ == "__main__":
    import asyncio

    loop = asyncio.new_event_loop()
    loop.run_until_complete(main_inference())
