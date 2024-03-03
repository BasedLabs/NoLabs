import json
from typing import List, Tuple

from nolabs.api_models.biobuddy import SendMessageRequest, SendMessageResponse
from nolabs.api_models.biobuddy import Message as ApiMessage
from nolabs.domain.experiment import ExperimentId
from nolabs.features.biobuddy.data_models.message import Message
from nolabs.features.biobuddy.file_management import FileManagement

import biobuddy_microservice
from biobuddy_microservice import DefaultApi
from biobuddy_microservice import ApiClient
from biobuddy_microservice import Configuration

import rcsb_pdb_query_microservice
from rcsb_pdb_query_microservice import DefaultApi as rcsbApiDefaultApi
from rcsb_pdb_query_microservice import ApiClient as rcsbApiClient
from rcsb_pdb_query_microservice import Configuration as rcsbApiConfiguration

import pubmed_query_microservice
from pubmed_query_microservice import DefaultApi as pubmedApiDefaultApi
from pubmed_query_microservice import ApiClient as pubmedApiClient
from pubmed_query_microservice import Configuration as pubmedApiConfiguration

import chembl_query_microservice
from chembl_query_microservice import DefaultApi as chemblApiDefaultApi
from chembl_query_microservice import ApiClient as chemblApiClient
from chembl_query_microservice import Configuration as chemblApiConfiguration

from nolabs.features.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.infrastructure.settings import Settings
from nolabs.utils.fasta import create_upload_file_from_string

tools = [
    {
        "type": "function",
        "function": {
            "name": "query_rcsb_pdb_by_id",
            "description": "Query RCSB PDB by protein ids",
            "parameters": {
                "type": "object",
                "properties": {
                    "pdb_ids": {
                        "type": "array",
                        "description": "Query rcsb pdb by protein ids",
                        "items": {
                            "type": "string"
                        }
                    },
                },
                "required": ["pdb_ids"],
            },
        }
    },
    {
        "type": "function",
        "function": {
            "name": "query_pubmed",
            "description": "Query Pubmed database by search string",
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "The string of the search query which is used to search pubmed"
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "The maximum number of articles to return"
                    },
                },
                "required": ["query", "max_results"],
            },
        }
    },
    {
            "type": "function",
            "function": {
                "name": "query_chembl",
                "description": "Query Chembl database by search term to find the molecule",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "search_term": {
                            "type": "string",
                            "description": "The search term of the molecule for which we want to find the results"
                        },
                    },
                    "required": ["search_term"],
                },
            }
        }
]

class SendMessageFeature:
    def __init__(self, settings: Settings, file_management: FileManagement,
                 targets_file_management: TargetsFileManagement):
        self._file_management = file_management
        self._targets_file_management = targets_file_management
        self._settings = settings

    def handle(self, request: SendMessageRequest) -> SendMessageResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id)
        messages = self._file_management.load_conversation(experiment_id)

        configuration = Configuration(
            host=self._settings.biobuddy_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            request = biobuddy_microservice.SendMessageToBioBuddyRequest(message_content=request.message_content,
                                                                         previous_messages=[{'role': message.role,
                                                                                             'content': message.content}
                                                                                            for message in messages],
                                                                         tools=tools)
            assistant_message = api_instance.predict_send_message_post(send_message_to_bio_buddy_request=request).chatgpt_reply

            if assistant_message.tool_calls and assistant_message.tool_calls.actual_instance:
                if assistant_message.tool_calls.actual_instance[0]['function']['name'] == "query_rcsb_pdb_by_id":
                    ids = json.loads(assistant_message.tool_calls.actual_instance[0]['function']['arguments'])["pdb_ids"]
                    results = self.query_rcsb_pdb_by_id(ids)
                    reply_content = f"Biobuddy called query_rcsb_pdb_by_id(), result: \n {results}"

                    for result in results:
                        temp_res = result[0][:]
                        file = create_upload_file_from_string(temp_res, "protein.fasta")
                        self._targets_file_management.store_target(experiment_id, file)

                    self._file_management.update_conversation(experiment_id,
                                                              Message(role='user', content=request.message_content,
                                                                      type='text'))
                    self._file_management.update_conversation(experiment_id,
                                                              Message(role='assistant',
                                                                      content=str(reply_content),
                                                                      type='text'))
                elif assistant_message.tool_calls.actual_instance[0]['function']['name'] == "query_pubmed":
                    query = json.loads(assistant_message.tool_calls.actual_instance[0]['function']['arguments'])[
                        "query"]
                    max_results = json.loads(assistant_message.tool_calls.actual_instance[0]['function']['arguments'])[
                        "max_results"]

                    articles = self.query_pubmed(query, max_results)

                    print(articles)


            else:
                self._file_management.update_conversation(experiment_id,
                                                          Message(role='user', content=request.message_content.actual_instance,
                                                                  type='text'))
                self._file_management.update_conversation(experiment_id,
                                                          Message(role='assistant',
                                                                  content=str(assistant_message.content.actual_instance),
                                                                  type='text'))


        return SendMessageResponse(biobuddy_response=ApiMessage(role='user', content=str(assistant_message.content), type='text'))

    def query_rcsb_pdb_by_id(self, pdb_ids: List[str]) -> List[Tuple[str, str]]:
        configuration = rcsbApiConfiguration(
            host=self._settings.rcsb_pdb_query_host,
        )
        with rcsbApiClient(configuration=configuration) as client:
            api_instance = rcsbApiDefaultApi(client)
            request = rcsb_pdb_query_microservice.GetFastaFilesByIdsRequest(rcsb_pdb_ids=pdb_ids)
            response = api_instance.predict_fetch_fastas_by_ids_post(get_fasta_files_by_ids_request=request)

            return [(result.fasta_contents, result.link) for result in response.fasta_contents]

    def query_pubmed(self, query: str, max_results: int) -> List[pubmed_query_microservice.FetchedArticle]:
        configuration = pubmedApiConfiguration(
            host=self._settings.pubmed_query_host,
        )
        with pubmedApiClient(configuration=configuration) as client:
            api_instance = pubmedApiDefaultApi(client)
            request = pubmed_query_microservice.PubMedSearchRequest(search_terms=query, max_results=max_results)
            response = api_instance.search_search_pubmed_articles_post(pub_med_search_request=request)

            return response.articles

    def query_chembl(self, search_term: str) -> List[chembl_query_microservice.Molecule]:
        configuration = chemblApiConfiguration(
            host=self._settings.chembl_query_host,
        )
        with chemblApiClient(configuration=configuration) as client:
            api_instance = chemblApiDefaultApi(client)
            request = chembl_query_microservice.ChEMBLMoleculeRequest(search_term=search_term)
            response = api_instance.query_query_chembl_post(ch_embl_molecule_request=request)

            return response.molecules


