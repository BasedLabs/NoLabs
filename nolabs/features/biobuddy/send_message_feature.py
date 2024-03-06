import json
from typing import List, Tuple

from nolabs.api_models.biobuddy import SendMessageRequest, SendMessageResponse
from nolabs.api_models.biobuddy import Message as ApiMessage
from nolabs.domain.experiment import ExperimentId
from nolabs.features.biobuddy.data_models.message import Message, FunctionCall, FunctionParam, RegularMessage
from nolabs.api_models.biobuddy import RegularMessage as ApiRegularMessage
from nolabs.api_models.biobuddy import FunctionCall as ApiFunctionCall
from nolabs.api_models.biobuddy import FunctionParam as ApiFunctionParam
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

from nolabs.features.drug_discovery.services.ligand_file_management import LigandsFileManagement
from nolabs.features.drug_discovery.services.target_file_management import TargetsFileManagement
from nolabs.infrastructure.settings import Settings
from nolabs.utils.fasta import create_upload_file_from_string
from nolabs.utils.sdf import smiles_to_sdf_string

tools = [
    {
        "type": "function",
        "function": {
            "name": "query_rcsb_pdb_by_id",
            "description": "Query RCSB PDB by protein ids. Don't use DNA or RNA ids. If a user asks for to pull targets/proteins then invoke this method.",
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
            "name": "query_rcsb_pdb_by_protein_names",
            "description": "Query RCSB PDB by protein  names. "
                           "If a user asks for to pull targets/proteins but does not specify ids then invoke this "
                           "method. Wathch out for the plural nouns, i.e. if a user asks to pull rhodopsins, then query rhodopsin.",
            "parameters": {
                "type": "object",
                "properties": {
                    "protein_names": {
                        "type": "array",
                        "description": "Query rcsb pdb by protein names",
                        "items": {
                            "type": "string"
                        }
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "Number of proteins to pull. If a user asks to pull exactly all of the "
                                       "proteins from a certain query,"
                                       "then don't add this parameter. If the user asks to pull some "
                                       "targets/proteins, then set this parameter to 10 by default.",
                    },
                },
                "required": ["protein_names"],
            },
        }
    },
    {
        "type": "function",
        "function": {
            "name": "query_pubmed",
            "description": "Query Pubmed database by search string. Query only if the user asks for articles/latest research papers.",
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
                 targets_file_management: TargetsFileManagement,
                 ligand_file_management: LigandsFileManagement):
        self._file_management = file_management
        self._targets_file_management = targets_file_management
        self._ligand_file_management = ligand_file_management
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

            previous_messages = []
            for message in messages:
                if message.type == "text":
                    previous_messages.append({'role': message.role,
                                          'content': str(message.message.content)})
                elif message.type == "function":
                    previous_messages.append({'role': message.role,
                                              'content': f"I called {message.message.function_name} with parameters:"
                                                         f"{message.message.parameters}"})


            request = biobuddy_microservice.SendMessageToBioBuddyRequest(message_content=request.message_content,
                                                                         previous_messages=[{'role': message.role,
                                                                                             'content': str(message.message)}
                                                                                            for message in messages],
                                                                         tools=tools)
            assistant_message = api_instance.predict_send_message_post(
                send_message_to_bio_buddy_request=request).chatgpt_reply

            print(assistant_message)
            if assistant_message.tool_calls and assistant_message.tool_calls.actual_instance:
                function_name = assistant_message.tool_calls.actual_instance[0]['function']['name']
                arguments = json.loads(assistant_message.tool_calls.actual_instance[0]['function']['arguments'])
                if function_name == "query_rcsb_pdb_by_id":
                    pdb_ids = arguments["pdb_ids"]
                    results = self.query_rcsb_pdb_by_id(pdb_ids)
                    reply_content = f"Sure! Pulled some targets, check them out: \n"
                    for result in results:
                        reply_content += f'{result[1]} \n'
                    for result in results:
                        temp_res = result[0][:]
                        file = create_upload_file_from_string(temp_res, "protein.fasta")
                        additional_metadata = {
                            "link": result[1]
                        }
                        self._targets_file_management.store_target(experiment_id, file, additional_metadata)

                    self._file_management.update_conversation(experiment_id,
                                                              Message(role='user',
                                                                      message=RegularMessage(
                                                                          content=request.message_content
                                                                      ),
                                                                      type='text')
                                                              )
                    self._file_management.update_conversation(experiment_id,
                                                              Message(role='assistant',
                                                                      message=FunctionCall(
                                                                          function_name="query_rcsb_pdb_by_id",
                                                                          parameters=[
                                                                              FunctionParam(name="pdb_ids",
                                                                                            value=pdb_ids)
                                                                          ]),
                                                                      type="function")
                                                              )
                    return SendMessageResponse(biobuddy_response=ApiMessage(role='assistant',
                                                                            message=ApiFunctionCall(
                                                                                function_name="query_rcsb_pdb_by_id",
                                                                                parameters=[
                                                                                    ApiFunctionParam(name="pdb_ids",
                                                                                                  value=pdb_ids)
                                                                                ]),
                                                                            type="function")
                                               )

                if function_name == "query_rcsb_pdb_by_protein_names":
                    protein_names = arguments["protein_names"]
                    max_results = None
                    if "max_results" in arguments:
                        max_results = arguments["max_results"]
                    results = self.query_rcsb_pdb_by_protein_names(protein_names, max_results)
                    reply_content = f"Sure! Pulled some targets, check them out: \n"
                    for result in results:
                        reply_content += f'{result[1]} \n'
                    for result in results:
                        temp_res = result[0][:]
                        file = create_upload_file_from_string(temp_res, "protein.fasta")
                        additional_metadata = {
                            "link": result[1]
                        }
                        self._targets_file_management.store_target(experiment_id, file, additional_metadata)

                    self._file_management.update_conversation(experiment_id,
                                                              Message(role='user',
                                                                      message=RegularMessage(
                                                                          content=request.message_content
                                                                      ),
                                                                      type='text')
                                                              )
                    self._file_management.update_conversation(experiment_id,
                                                              Message(role='assistant',
                                                                      message=FunctionCall(
                                                                          function_name="query_rcsb_pdb_by_protein_names",
                                                                          parameters=[
                                                                              FunctionParam(name="protein_names",
                                                                                            value=protein_names)
                                                                          ]),
                                                                      type="function")
                                                              )
                    return SendMessageResponse(biobuddy_response=ApiMessage(role='assistant',
                                                                            message=ApiFunctionCall(
                                                                                function_name="query_rcsb_pdb_by_protein_names",
                                                                                parameters=[
                                                                                    ApiFunctionParam(name="protein_names",
                                                                                                  value=protein_names)
                                                                                ]),
                                                                            type="function")
                                               )

                elif function_name == "query_pubmed":
                    query = arguments["query"]
                    max_results = arguments["max_results"]
                    articles = self.query_pubmed(query, max_results)
                    reply_content = f"Sure! Pulled some ligands, check them out: \n"
                    for article in articles:
                        reply_content += f"{article.title} \n"

                    self._file_management.update_conversation(experiment_id,
                                                              Message(role='user',
                                                                      message=RegularMessage(
                                                                          content=request.message_content
                                                                      ),
                                                                      type='text'))
                    self._file_management.update_conversation(experiment_id,
                                                              Message(role='assistant',
                                                                      message=FunctionCall(
                                                                          function_name="query_pubmed",
                                                                          parameters=[
                                                                              FunctionParam(name="query",
                                                                                            value=query),
                                                                              FunctionParam(name="max_results",
                                                                                            value=max_results),
                                                                          ]),
                                                                      type="function")
                                                              )
                    return SendMessageResponse(biobuddy_response=ApiMessage(role='assistant',
                                                                            message=ApiFunctionCall(
                                                                                function_name="query_pubmed",
                                                                                parameters=[
                                                                              ApiFunctionParam(name="query",
                                                                                            value=query),
                                                                              ApiFunctionParam(name="max_results",
                                                                                            value=max_results),
                                                                          ]),
                                                                      type="function"))

                elif function_name == "query_chembl":
                    search_term = arguments["search_term"]
                    molecules = self.query_chembl(search_term)
                    reply_content = f"Sure! Pulled some ligands, check them out: \n"
                    for molecule in molecules:
                        reply_content += f'{molecule.chembl_id} \n'
                        molecule_smiles = molecule.smiles
                        sdf_content = smiles_to_sdf_string(molecule_smiles)
                        file = create_upload_file_from_string(sdf_content, f"{molecule.chembl_id}.sdf")
                        additional_metadata = {
                            "link": molecule.link
                        }
                        self._ligand_file_management.store_lone_ligand(experiment_id, file, additional_metadata)
                    self._file_management.update_conversation(experiment_id,
                                                              Message(role='user',
                                                                      message=RegularMessage(
                                                                          content=request.message_content
                                                                      ),
                                                                      type='text')
                                                              )
                    self._file_management.update_conversation(experiment_id,
                                                              Message(role='assistant',
                                                                      message=FunctionCall(
                                                                          function_name='query_chembl',
                                                                          parameters=[
                                                                              FunctionParam(name="search_term",
                                                                                            value=search_term),
                                                                          ]
                                                                      ),
                                                                      type="function"
                                                                      )
                                                              )
                    return SendMessageResponse(biobuddy_response=ApiMessage(role='assistant',
                                                                            message=ApiFunctionCall(
                                                                                function_name='query_chembl',
                                                                                parameters=[
                                                                                    ApiFunctionParam(name="search_term",
                                                                                                     value=search_term),
                                                                                ]
                                                                            ),
                                                                            type='function')
                                               )
            else:
                self._file_management.update_conversation(experiment_id,
                                                          Message(role='user',
                                                                  message=RegularMessage(
                                                                      content=request.message_content
                                                                  ),
                                                                  type='text')
                                                          )
                self._file_management.update_conversation(experiment_id,
                                                          Message(role='assistant',
                                                                  message=RegularMessage(
                                                                      content=str(
                                                                          assistant_message.content.actual_instance)
                                                                  ),
                                                                  type='text')
                                                          )
                return SendMessageResponse(biobuddy_response=ApiMessage(role='assistant',
                                                                        message=ApiRegularMessage(
                                                                            content=str(assistant_message.content.actual_instance)
                                                                        ),
                                                                        type='text')
                                           )



    def query_rcsb_pdb_by_id(self, pdb_ids: List[str]) -> List[Tuple[str, str]]:
        configuration = rcsbApiConfiguration(
            host=self._settings.rcsb_pdb_query_host,
        )
        with rcsbApiClient(configuration=configuration) as client:
            api_instance = rcsbApiDefaultApi(client)
            request = rcsb_pdb_query_microservice.GetFastaFilesByIdsRequest(rcsb_pdb_ids=pdb_ids)
            response = api_instance.fetch_fetch_fastas_by_ids_post(get_fasta_files_by_ids_request=request)

            return [(result.fasta_contents, result.link) for result in response.fasta_contents]

    def query_rcsb_pdb_by_protein_names(self, protein_names: List[str], max_results: int = None) -> List[Tuple[str, str]]:
        configuration = rcsbApiConfiguration(
            host=self._settings.rcsb_pdb_query_host,
        )
        results = []
        with rcsbApiClient(configuration=configuration) as client:
            api_instance = rcsbApiDefaultApi(client)
            for protein_name in protein_names:
                request = rcsb_pdb_query_microservice.GetFastaFilesBySearchQueryRequest(search_query=protein_name,
                                                                                        max_results=max_results)
                response = api_instance.fetch_fetch_fastas_by_search_query_post(get_fasta_files_by_search_query_request=request)
                for result in response.fasta_contents:
                    results.append((result.fasta_contents, result.link))
        return results

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
