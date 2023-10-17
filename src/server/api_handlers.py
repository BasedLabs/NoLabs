import time

from flask import Request
import src.server.services.inference_service as inference_service
from src.server import settings
from src.server.services.fasta_reader import get_sequences
from src.server.services.oboreader import read_obo


class ApiHandler:
    pass


class AminoAcidLabApiHandler(ApiHandler):
    def inference(self, request: Request) -> dict:
        amino_acid_input_sequence = request.form['inputSequence']
        amino_acid_input_sequence_files = request.files.getlist('inputSequenceFile')
        if not amino_acid_input_sequence:
            amino_acid_input_sequence = \
                [seq for seq in get_sequences(amino_acid_input_sequence_files)][0]
        pipeline = inference_service.create_pipeline(use_gpu=True, is_test=settings.is_test)
        localisation_result = inference_service.get_localisation_output(pipeline=pipeline,
                                                                        amino_acid_sequence=amino_acid_input_sequence)

        folding_result = inference_service.get_folding_output(pipeline=pipeline,
                                                              amino_acid_sequence=amino_acid_input_sequence)
        gene_ontology_result = inference_service.get_gene_ontology_output(pipeline=pipeline,
                                                                          amino_acid_sequence=amino_acid_input_sequence)
        solubility = inference_service.get_solubility_output(pipeline=pipeline,
                                                             amino_acid_sequence=amino_acid_input_sequence)

        esm_protein_localization = {key: value for key, value in localisation_result}
        obo_graph = read_obo(gene_ontology_result)
        return {
            'sequence': amino_acid_input_sequence,
            'localisation': {
                'mithochondria': esm_protein_localization["Mitochondrial Proteins"],
                'nucleus': esm_protein_localization["Nuclear Proteins"],
                'cytoplasm': esm_protein_localization["Cytosolic Proteins"],
                'other': esm_protein_localization["Other Proteins"],
                'extracellular': esm_protein_localization["Extracellular/Secreted Proteins"]
            },
            'folding': folding_result,
            'oboGraph': obo_graph,
            'solubility': solubility
        }


class AminoAcidLabApiMockHandler(AminoAcidLabApiHandler):
    def inference(self, request: Request) -> dict:
        amino_acid_input_sequence = request.form['inputSequence']
        amino_acid_input_sequence_files = request.files.getlist('inputSequenceFile')

        return {
            'sequence': 'AAAAAAAAAA',
            'localisation': {
                'mithochondria': 0.2,
                'nucleus': 0.5,
                'cytoplasm': 0.1,
                'other': 0.4,
                'extracellular': 0.3,
            },
            'folding': open('mock_data/test.pdb').read(),
            'oboGraph': {
                'GO:123': {'name': 'GO:123', 'namespace': 'biological_process', 'edges': {}},
                'GO:234': {'name': 'GO:234', 'namespace': 'biological_process', 'edges': {}}
            },
            'solubility': 0.5
        }


class DrugTargetApiHandler(ApiHandler):
    def inference(self, request):
        ligand_files = request.files.getlist('sdfFileInput')
        protein_files = request.files.getlist('proteinFileInput')

        pipeline = inference_service.create_pipeline(use_gpu=settings.use_gpu, is_test=settings.is_test)

        inference_service.save_uploaded_files(pipeline, ligand_files)
        inference_service.save_uploaded_files(pipeline, protein_files)

        inference_service.generate_dti_results(pipeline, ligand_files, protein_files)
        pdb_content, protein_name, ligands_sdf_contents, ligand_names, affinity_list \
            = inference_service.get_dti_results(pipeline, ligand_files)

        if not ligand_files:
           return 'You must provide either smiles text or smiles files', 400

        if not protein_files:
           return 'You must provide a protein .pdb file', 400

        res = {'name': request.form['experimentName'], 'data': [{'proteinName': protein_name,
               'ligandName': ligand_name,
               'pdb': pdb_content,
               'sdf': ligand_content,
               'affinity': affinity} for ligand_name, ligand_content, affinity \
               in zip(ligand_names, ligands_sdf_contents, affinity_list)]}

        return res

    def get_experiments(self):
        return []

    def get_experiment(self, request):
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_name = request.args.get('name')

        return {'name': experiment_name, 'data': []}

    def delete_experiment(self):
        # Get name of experiment here and delete it based on this name
        return 200


class DrugTargetApiMockHandler(DrugTargetApiHandler):
    def inference(self, request):
        ligand_files = request.files.getlist('sdfFileInput')
        protein_files = request.files.getlist('proteinFileInput')
        experiment_name = request.form['experimentName']

        return {'name': experiment_name, 'data': [{
            'proteinName': "AHAHAHAHAHHAHA2222222222222222",
            'ligandName': 'LALSDLASDLASLDASLDA IAM CRAZYYYY',
            'pdb': open('mock_data/test.pdb').read(),
            'sdf': open('mock_data/test.sdf').read(),
            'affinity': 10
        }]}

    def get_experiments(self):
        experiments = [
            {'name': 'Experiment 10'},
            {'name': 'Experiment 11'}
        ]
        return experiments

    def get_experiment(self, request):
        time.sleep(10)
        # get name of the experiment and get EXISTING SAVED data based on this name
        experiment_name = request.args.get('name')

        return {'name': experiment_name, 'data': [{
            'proteinName': "AHAHAHAHAHHAHA",
            'ligandName': 'LALSDLASDLASLDASLDA IAM CRAZYYYY',
            'pdb': open('mock_data/test.pdb').read(),
            'sdf': open('mock_data/test.sdf').read(),
            'affinity': 10
        }]}

    def delete_experiment(self):
        # Get name of experiment here and delete it based on this name
        return 200
