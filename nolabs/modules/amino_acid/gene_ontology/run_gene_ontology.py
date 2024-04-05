import numbers
import os
import pprint

# noinspection PyUnresolvedReferences
import obonet
import pickle
from typing import Dict, List, Tuple

from gene_ontology_microservice import Configuration, ApiClient, DefaultApi, RunGeneOntologyPredictionRequest

from nolabs.api_models.amino_acid.common_models import RunAminoAcidRequest
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.amino_acid import AminoAcid
from nolabs.domain.experiment import ExperimentId, ExperimentName
from nolabs.modules.amino_acid.gene_ontology.services.file_management import FileManagement
from nolabs.modules.amino_acid.run_aa_inference_feature_base import RunAminoAcidInferenceFeature
from nolabs.infrastructure.settings import Settings
from nolabs.domain.gene_ontology import OboNode
from nolabs.api_models.amino_acid.gene_ontology import RunGeneOntologyResponse, AminoAcidResponse, \
    RunGeneOntologyResponseDataNode
from nolabs.utils import generate_uuid
from nolabs.utils.fasta import FastaReader


class RunGeneOntologyFeature(RunAminoAcidInferenceFeature[FileManagement]):
    def __init__(self, settings: Settings, file_management: FileManagement):
        super().__init__(file_management)
        self._settings = settings

    def _read_obo(self, obo: List[Tuple[str, float]]) -> Dict[str, OboNode]:
        ids_dict = {name: prob for name, prob in obo}
        ids_dict = {key: value for key, value in ids_dict.items() if value >= 0.7}
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        pickled_file = os.path.join(__location__, 'services', 'go.pickle')

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

    async def handle(self,
                     request: RunAminoAcidRequest) -> RunGeneOntologyResponse:
        assert request

        experiment_id = ExperimentId(request.experiment_id) if request.experiment_id else ExperimentId(generate_uuid())

        await self._setup_experiment(experiment_id, request)

        configuration = Configuration(
            host=self._settings.gene_ontology_host,
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            amino_acids: List[AminoAcid] = []

            if request.amino_acid_sequence:
                amino_acids.append(AminoAcid(name='unknown', sequence=request.amino_acid_sequence))

            if request.fastas:
                fasta_reader = FastaReader()
                for fasta in request.fastas:
                    fasta_content = await fasta.read()
                    for amino_acid in fasta_reader.get_ids2seqs(fasta_content.decode('utf-8')):
                        if amino_acid.name not in [aa.name for aa in amino_acids]:
                            amino_acids.append(amino_acid)

            results: List[AminoAcidResponse] = []
            if not amino_acids:
                raise NoLabsException(['No amino acids'], ErrorCodes.no_amino_acids)
            for amino_acid in amino_acids:
                result = api_instance.run_go_prediction_run_go_prediction_post(
                    run_gene_ontology_prediction_request=RunGeneOntologyPredictionRequest(
                        amino_acid_sequence=amino_acid.sequence
                    )
                )

                if result.errors:
                    raise NoLabsException(result.errors, ErrorCodes.amino_acid_solubility_run_error)

                confidences = [(g.name, g.confidence) for g in result.go_confidence]

                results.append(AminoAcidResponse(
                    sequence=amino_acid.sequence,
                    name=amino_acid.name,
                    go={key: RunGeneOntologyResponseDataNode(
                        name=value.name,
                        namespace=value.namespace,
                        edges=value.edges
                    ) for key, value in self._read_obo(confidences).items()}
                ))

            self._file_management.set_result(experiment_id=experiment_id,
                                             data=results)

            metadata = self._file_management.get_metadata(experiment_id=experiment_id)
            return RunGeneOntologyResponse(experiment_id=experiment_id.value,
                                           experiment_name=metadata.name.value,
                                           amino_acids=results
                                           )
