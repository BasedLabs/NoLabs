from langchain.embeddings import OpenAIEmbeddings
from langchain.vectorstores import FAISS


class Memory:
    def __init__(self, **kwargs):
        self._embeddings = OpenAIEmbeddings()

    def get_embeddings(self):
        return self._embeddings
