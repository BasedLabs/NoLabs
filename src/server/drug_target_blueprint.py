from flask import Blueprint

import argparse

from flask import Flask, render_template, request
import src.server.services.inference_service as inference_service
from src.server.services.oboreader import read_obo
from src.server.services.fasta_reader import get_sequences

drug_target_bp = Blueprint('drug-target', __name__)


@drug_target_bp.route('/inference', methods=['POST'])
def inference():
    # ligand_files = request.files.getlist('sdfFileInput')
    # protein_files = request.files.getlist('proteinFileInput')

    # pipeline = inference_service.create_pipeline(use_gpu=use_gpu, is_test=is_test)

    # inference_service.save_uploaded_files(pipeline, ligand_files)
    # inference_service.save_uploaded_files(pipeline, protein_files)

    # inference_service.generate_dti_results(pipeline, ligand_files, protein_files)
    # pdb_content, protein_name, ligands_sdf_contents, ligand_names, affinity_list \
    #     = inference_service.get_dti_results(pipeline, ligand_files)

    # if not ligand_files:
    #    return 'You must provide either smiles text or smiles files', 400

    # if not protein_files:
    #    return 'You must provide a protein .pdb file', 400

    # res = [{'proteinName': protein_name,
    #        'ligandName': ligand_name,
    #        'pdb': pdb_content,
    #        'sdf': ligand_content,
    #        'affinity': affinity} for ligand_name, ligand_content, affinity \
    #        in zip(ligand_names, ligands_sdf_contents, affinity_list)]
    #
    # return res

    experiment_name = request.form['experimentName']

    # SAVE inference locally with experiment name

    return {'name': 'Experiment 10', 'data': [{
        'proteinName': "AHAHAHAHAHHAHA",
        'ligandName': 'LALSDLASDLASLDASLDA IAM CRAZYYYY',
        'pdb': open('test.pdb').read(),
        'sdf': open('test.sdf').read(),
        'affinity': 10
    }]}


@drug_target_bp.route('/experiments')
def get_experiments():
    experiments = [
        {'name': 'Experiment 10'},
        {'name': 'Experiment 11'}
    ]
    return experiments

@drug_target_bp.route('/load-experiment/<name>')
def get_experiment():
    # get name of the experiment and get EXISTING SAVED data based on this name

    return {'name': 'Experiment 10', 'data': [{
        'proteinName': "AHAHAHAHAHHAHA",
        'ligandName': 'LALSDLASDLASLDASLDA IAM CRAZYYYY',
        'pdb': open('test.pdb').read(),
        'sdf': open('test.sdf').read(),
        'affinity': 10
    }]}


@drug_target_bp.route('/delete-experiment', methods=['DELETE'])
def drug_target_discovery_delete_experiment():
    # Get name of experiment here and delete it based on this name
    return 200
