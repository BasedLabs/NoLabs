# Arxiv AI abstracts search

This microservice contains LLM RAG search over arxiv abstracts.
How to use this docker image

## Getting started

1) Generate a new token for docker registry
https://github.com/settings/tokens/new?scopes=read:packages
Select 'read:packages' (should be automatically selected when navigating link above).
2) Download ChromaDb
```bash
$ cd ../..
$ make download-arxiv-abstracts-db 
```
3) Start docker
```bash
$ docker login ghcr.io -u username -p ghp_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
$ docker compose up arxiv-ai-abstractions-search
```

If you're running api.py then you can access UI in browser `http://0.0.0.0:8001/chat`