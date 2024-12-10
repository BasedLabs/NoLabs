import json
import logging
import os
import uuid
from threading import Lock

from dotenv import load_dotenv

load_dotenv(".env")

from functools import lru_cache
from typing import List, AsyncGenerator, Iterable, Sequence, Union
from settings import settings

from api_models import ChatResponse, DocumentContext, MessageType, ChatRequest
from langchain_core.documents import Document
from langgraph.graph import MessagesState
from logging import getLogger
from chromadb.api.types import IncludeEnum
from langchain_huggingface import HuggingFaceEmbeddings
from langchain_chroma import Chroma
from paperqa import EmbeddingModel, Docs, Settings
from paperqa.llms import VectorStore
from paperqa.types import Embeddable, Doc, Text

HISTORY_FILE = "history.json"
file_lock = Lock()  # Lock to ensure thread safety


@lru_cache
def get_vector_store():
    embeddings = HuggingFaceEmbeddings(model_name="sentence-transformers/all-mpnet-base-v2")
    return (Chroma(
        collection_name="arxiv_abstracts",
        persist_directory=settings.chroma_db_path,  # This should be the same directory as before
        embedding_function=embeddings
    ), embeddings)


class ChromaDbVector(VectorStore):
    texts_hashes: set[int] = [0]

    def __contains__(self, item) -> bool:
        return True

    def __len__(self) -> int:
        return len(self.texts_hashes)

    def clear(self) -> None:
        pass

    def add_texts_and_embeddings(self, texts: Iterable[Embeddable]) -> None:
        pass

    async def similarity_search(
            self, query: str, k: int, embedding_model: EmbeddingModel
    ) -> tuple[Sequence[Embeddable], list[float]]:
        chroma, embeddings = get_vector_store()
        e = embeddings.embed_documents([query])
        res = chroma._collection.query(e, include=[IncludeEnum.embeddings, IncludeEnum.metadatas, IncludeEnum.documents,
                                                   IncludeEnum.distances], n_results=k)
        texts = []

        for emb, met, d in zip(res['embeddings'][0], res['metadatas'][0], res['documents'][0]):
            doc = Doc(embedding=emb, docname=met['title'], citation=met['id'], dockey=met['id'])
            text = Text(embeddings=emb, text=d, name=met['title'], doc=doc)
            texts.append(text)

        return (
            texts,
            res['distances'][0]
        )


@lru_cache
def get_docs():
    return Docs(texts_index=ChromaDbVector())


@lru_cache
def get_settings():
    settings = Settings(embedding='st-all-mpnet-base-v2')
    settings.answer.answer_max_sources = 3
    settings.answer.evidence_k = 5
    settings.prompts.use_json = False
    return settings


print('Initializing vector store, it can take a minute or two')

get_vector_store()
get_docs()


def save_history(req_res: Union[ChatRequest, ChatResponse]):
    # Ensure the file exists
    if not os.path.exists(HISTORY_FILE):
        with open(HISTORY_FILE, "w") as f:
            json.dump([], f)  # Initialize an empty list

    with file_lock:
        with open(HISTORY_FILE, "r") as f:
            history = json.load(f)

        history.append(req_res.model_dump())

        with open(HISTORY_FILE, "w") as f:
            json.dump(history, f, indent=4)


async def inference(request: ChatRequest) -> AsyncGenerator[ChatResponse, None]:
    save_history(request)

    docs = get_docs()
    session = await docs.aquery(
        query=request.message,
        settings=get_settings(),
    )

    documents = []

    for context in session.contexts:
        document = DocumentContext(
            page_content=context.text.text,
            title=context.text.name,
            id=context.text.doc.dockey
        )
        documents.append(document)

    response = ChatResponse(
        id=str(session.id),
        message_type=MessageType.ai,
        context=documents,
        content=session.answer
    )
    save_history(response)
    yield response


def set_openai_api_key(api_key: str):
    settings.openai_api_key = api_key
    os.environ['OPENAI_API_KEY'] = api_key


async def history() -> List[ChatResponse]:
    if not os.path.exists(HISTORY_FILE):
        return []

    j = json.load(open(HISTORY_FILE, "r"))
    result = []
    for item in j:
        if 'message' in item:
            result.append(
                ChatResponse(id=str(uuid.uuid4()), context=[], message_type=MessageType.human, content=item['message']))
        else:
            result.append(ChatResponse(**item))
    return result


def clear_history():
    os.remove(HISTORY_FILE)


if __name__ == '__main__':
    import asyncio

    loop = asyncio.new_event_loop()
    loop.run_until_complete(inference(request=ChatRequest(message='Tell me about quarks flavours')))
