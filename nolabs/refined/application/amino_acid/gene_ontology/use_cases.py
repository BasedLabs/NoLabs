__all__ = [
    'GetGeneOntologyJobFeature',
    'RunGeneOntologyFeature'
]

import numbers
import os
import pickle
import uuid
from typing import List, Dict, Any, Tuple
from uuid import UUID

from gene_ontology_microservice import DefaultApi, RunGeneOntologyPredictionRequest

from nolabs.domain.gene_ontology import OboNode
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.amino_acid.api_models import RunAminoAcidRequest
from nolabs.refined.application.amino_acid.gene_ontology.api_models import GetJobResponse, \
    JobPropertiesResponse, \
    AminoAcidResponse, JobFastaPropertyResponse, RunJobResponse, RunGeneOntologyResponseDataNode
from nolabs.refined.application.amino_acid.services import get_input_proteins
from nolabs.refined.domain.models.common import JobId, Experiment, ExperimentId, \
    JobName, Protein
from nolabs.refined.domain.models.gene_ontology import GeneOntologyJob
from nolabs.refined.infrastructure.settings import Settings


def map_to_amino_acid_response(protein: Protein) -> AminoAcidResponse:
    return AminoAcidResponse(
        sequence=protein.get_fasta(),
        name=str(protein.name),
        go=protein.gene_ontology
    )


class GetGeneOntologyJobFeature:
    async def handle(self, job_id: UUID) -> GetJobResponse:
        job_id = JobId(job_id)
        job: GeneOntologyJob = GeneOntologyJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        proteins = job.proteins

        return GetJobResponse(
            job_id=job_id.value,
            job_name=str(job.name),
            amino_acids=[
                map_to_amino_acid_response(protein)
                for protein in job.proteins
            ],
            properties=JobPropertiesResponse(
                fastas=[
                    JobFastaPropertyResponse(
                        filename=f.name.fasta_name,
                        content=f.get_fasta()
                    ) for f in proteins
                ]
            )
        )


class RunGeneOntologyFeature:
    _settings: Settings
    _api: DefaultApi

    def __init__(self, api: DefaultApi, settings: Settings):
        self._settings = settings
        self._api = api

    async def handle(self, request: RunAminoAcidRequest) -> RunJobResponse:
        assert request
        assert request.experiment_id

        if not request.fastas:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        experiment_id = ExperimentId(request.experiment_id)
        experiment: Experiment = Experiment.objects.with_id(experiment_id.value)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        if not request.job_id:
            job = GeneOntologyJob(
                id=JobId(uuid.uuid4()),
                name=JobName('New job'),
                experiment=experiment
            )
            job.save()
            job_id = JobId(job.id)
        else:
            job_id = JobId(request.job_id)
            job: GeneOntologyJob = GeneOntologyJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        try:
            input_proteins = await get_input_proteins(experiment=experiment, request=request)

            job.set_proteins(input_proteins)
            job.save(cascade=True)

            job.clear_result()

            results: List[AminoAcidResponse] = []

            for protein in job.proteins:
                result = self._api.run_go_prediction_run_go_prediction_post(
                    run_gene_ontology_prediction_request=RunGeneOntologyPredictionRequest(
                        amino_acid_sequence=protein.get_fasta()
                    )
                )

                if result.errors:
                    raise NoLabsException(ErrorCodes.gene_ontology_run_error, result.errors)

                confidences = [(g.name, g.confidence) for g in result.go_confidence]

                go = {key: RunGeneOntologyResponseDataNode(
                        name=value.name,
                        namespace=value.namespace,
                        edges=value.edges
                    ) for key, value in self._read_obo(confidences).items()}

                results.append(AminoAcidResponse(
                    sequence=protein.get_fasta(),
                    name=protein.name.value,
                    go=go
                ))

                job.set_result(protein, go)
                protein.set_gene_ontology(gene_ontology=go)
                protein.save()

            job.save(cascade=True)

            results: List[AminoAcidResponse] = []

            for protein in job.proteins:
                results.append(map_to_amino_acid_response(protein))

            return RunJobResponse(job_id=job_id.value,
                                  amino_acids=results)
        except Exception as e:
            print(e)
            if not isinstance(e, NoLabsException):
                raise NoLabsException(ErrorCodes.unknown_gene_ontology_error) from e
            raise e

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
