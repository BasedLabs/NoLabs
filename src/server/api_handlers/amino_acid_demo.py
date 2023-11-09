from flask import Request

from src.server.services.oboreader import read_obo
from src.server.services.experiments_structure_loader import ProteinLabExperimentsLoader, DTILabExperimentsLoader
from src.server.api_handlers.api_handler import ApiHandler


class AminoAcidLabDemoApiHandler(ApiHandler):
    def __init__(self):
        self.experiments_loader = ProteinLabExperimentsLoader()

    def get_experiments(self):
        return self.experiments_loader.load_experiments()

    def get_experiment(self, request):
        experiment_id = request.args.get('id')
        experiment_name = request.args.get('name')

        result = self.experiments_loader.load_experiment(experiment_id)
        localisation_result = result['localisation']
        folding_result = result['folding']
        gene_ontology_result = result['gene_ontology']
        solubility_result = result['solubility']

        obo_graph = read_obo(gene_ontology_result)

        return {'id': experiment_id, 'name': experiment_name, 'data': {
            'localisation': {
                'mithochondria': localisation_result["Mitochondrial Proteins"],
                'nucleus': localisation_result["Nuclear Proteins"],
                'cytoplasm': localisation_result["Cytosolic Proteins"],
                'other': localisation_result["Other Proteins"],
                'extracellular': localisation_result["Extracellular/Secreted Proteins"]
            },
            'folding': folding_result,
            'oboGraph': obo_graph,
            'solubility': solubility_result['solubility']
        }}

    def change_experiment_name(self, request: Request):
        return {'result': 'Dont look here! It is a demo'}

    def delete_experiment(self, request: Request):
        return {'result': 'Please no'}

    def inference(self, request: Request) -> dict:
        return {'result': 'This is a demo, deploy docker container or buy us a fancy gpu'}
