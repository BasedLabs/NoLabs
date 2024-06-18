__all__ = [
    'GetJobFeature',
    'RunJobFeature',
    'SetupJobFeature'
]

import numbers
import os
import pickle
from typing import List, Dict, Any, Tuple
from uuid import UUID

from gene_ontology_microservice import DefaultApi, RunGeneOntologyPredictionRequest
from mongoengine import Q

from nolabs.domain.gene_ontology import OboNode
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.application.use_cases.gene_ontology.api_models import RunGeneOntologyResponseDataNode, JobResponse, \
    JobResult, SetupJobRequest
from nolabs.domain.models.common import JobId, Experiment, JobName, Protein
from nolabs.domain.models.gene_ontology import GeneOntologyJob
from nolabs.utils import generate_uuid


def map_job_to_response(job: GeneOntologyJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        protein_ids=[p.iid.value for p in job.proteins],
        result=[
            JobResult(
                protein_id=item.protein_id,
                go={key: RunGeneOntologyResponseDataNode(name=obo['name'], namespace=obo['namespace'],
                                                         edges=obo['edges']) for key, obo in item.gene_ontology.items()}
            )
            for item in job.gene_ontologies
        ],
        experiment_id=job.experiment.id
    )


class GetJobFeature:
    """
    Use case - get job information.
    """

    async def handle(self, job_id: UUID) -> JobResponse:
        job_id = JobId(job_id)
        job: GeneOntologyJob = GeneOntologyJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        return map_job_to_response(job)


class RunJobFeature:
    _api: DefaultApi

    def __init__(self, api: DefaultApi):
        self._api = api

    async def handle(self, job_id: UUID) -> JobResponse:
        assert job_id

        job_id = JobId(job_id)
        job: GeneOntologyJob = GeneOntologyJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        result: List[Tuple[Protein, Dict[str, Any]]] = []

        for protein in job.proteins:
            response = self._api.run_go_prediction_run_go_prediction_post(
                run_gene_ontology_prediction_request=RunGeneOntologyPredictionRequest(
                    amino_acid_sequence=protein.get_amino_acid_sequence()
                )
            )

            if response.errors:
                raise NoLabsException(ErrorCodes.gene_ontology_run_error, response.errors)

            confidences = [(g.name, g.confidence) for g in response.go_confidence]

            go = {key: {'name': value.name, 'namespace': value.namespace, 'edges': value.edges} for key, value in
                  self._read_obo(confidences).items()}

            result.append((protein, go))

        job.set_result(result=result)
        job.save(cascade=True)

        for protein, go in result:
            protein.set_gene_ontology(go)
            protein.save()

        return map_job_to_response(job)


    def _read_obo(self, obo: List[Tuple[str, float]]) -> Dict[str, OboNode]:
        ids_dict = {name: prob for name, prob in obo}
        ids_dict = {key: value for key, value in ids_dict.items() if value >= 0.7}
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        pickled_file = os.path.join(__location__, 'go.pickle')

        # f = os.path.join(__location__, 'services', 'go.obo')
        # graph = obonet.read_obo(f)

        # with open(pickled_file, 'wb') as handle_w:
        #    pickle.dump(graph, handle_w, protocol=pickle.HIGHEST_PROTOCOL)
        with open(pickled_file, 'rb') as handle_r:
            graph = pickle.load(handle_r)

        ids = set(ids_dict.keys())

        def get_attributes(edge_data):
            for key in edge_data:
                if isinstance(key, numbers.Number):
                    return get_attributes(edge_data[key])
                return [k for k in edge_data.keys()]

        # cleanup a graph from the nodes that we don't have in list
        nodes_to_remove = []
        for node in graph.nodes:
            if node in ids:
                continue
            for in_edge in graph.in_edges(node):
                for out_edge in graph.out_edges(node):
                    if not graph.has_edge(in_edge[0], out_edge[1]):
                        attributes = graph.get_edge_data(out_edge[0], out_edge[1], 0)
                        if not attributes:
                            attributes = graph.get_edge_data(out_edge[0], out_edge[1])
                        graph.add_edge(in_edge[0], out_edge[1], **attributes)
            nodes_to_remove.append(node)

        for node in nodes_to_remove:
            graph.remove_node(node)

        shaped_graph: Dict[str, OboNode] = {}  # {nodeId: {name: '', edges: { nodeId: {linkType}, }}}
        for node in graph.nodes:
            if node not in shaped_graph:
                shaped_graph[node] = OboNode(name=graph.nodes[node]['name'],
                                             namespace=graph.nodes[node]['namespace'],
                                             edges={})
        for node in graph.nodes:
            for out_edge in graph.out_edges(node):
                out_node = out_edge[1]
                edge_data = get_attributes(graph.get_edge_data(node, out_edge[1]))
                shaped_graph[out_node].edges[node] = edge_data

        return shaped_graph


class SetupJobFeature:
    """
    Use case - create new or update existing job.
    """
    async def handle(self, request: SetupJobRequest) -> JobResponse:
        assert request

        job_id = JobId(request.job_id or generate_uuid())
        job_name = JobName(request.job_name or 'New gene ontology job')

        jobs: GeneOntologyJob = GeneOntologyJob.objects(Q(id=job_id.value) | Q(name=job_name.value))

        if not jobs:
            if not request.experiment_id:
                raise NoLabsException(ErrorCodes.invalid_experiment_id)

            experiment = Experiment.objects.with_id(request.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            job = GeneOntologyJob(
                id=job_id,
                name=job_name,
                experiment=experiment
            )
        else:
            job = jobs[0]

        proteins: List[Protein] = []
        for protein_id in request.proteins:
            protein = Protein.objects.with_id(protein_id)

            if not protein:
                raise NoLabsException(ErrorCodes.protein_not_found)

            proteins.append(protein)

        job.set_inputs(proteins=proteins)
        job.save(cascade=True)

        return map_job_to_response(job)
