from flask import Blueprint

import argparse

from flask import Flask, render_template, request
import src.server.services.inference_service as inference_service
from src.server.services.oboreader import read_obo
from src.server.services.fasta_reader import get_sequences

drug_target_bp = Blueprint('drug-target', __name__,
                        template_folder='templates')

@drug_target_bp.route('/')
def drug_target_discovery():
    return render_template('drug-target-page.html')


@drug_target_bp.route('/inference', methods=['POST'])
def inference_drug_target_discovery():
   #ligand_files = request.files.getlist('sdfFileInput')
   #protein_files = request.files.getlist('proteinFileInput')

   #pipeline = inference_service.create_pipeline(use_gpu=use_gpu, is_test=is_test)

   #inference_service.save_uploaded_files(pipeline, ligand_files)
   #inference_service.save_uploaded_files(pipeline, protein_files)


   #inference_service.generate_dti_results(pipeline, ligand_files, protein_files)
   #pdb_content, protein_name, ligands_sdf_contents, ligand_names, affinity_list \
   #     = inference_service.get_dti_results(pipeline, ligand_files)

   #if not ligand_files:
   #    return 'You must provide either smiles text or smiles files', 400

   #if not protein_files:
   #    return 'You must provide a protein .pdb file', 400

    #res = [{'proteinName': protein_name,
    #        'ligandName': ligand_name,
    #        'pdb': pdb_content,
    #        'sdf': ligand_content,
    #        'affinity': affinity} for ligand_name, ligand_content, affinity \
    #        in zip(ligand_names, ligands_sdf_contents, affinity_list)]
#
    #return {'drugTarget': res}

    return {'drugTarget': [{
        'proteinName': "AHAHAHAHAHHAHA",
        'ligandName': 'LALSDLASDLASLDASLDA IAM CRAZYYYY',
        'pdb': open('test.pdb').read(),
        'sdf': open('test.sdf').read(),
        'affinity': 10
    }]}


@drug_target_bp.route('/load-experiments')
def drug_target_load_experiments():
    experiments = [{
        'experimentName': 'Experiment 1',
        'proteinName': "AHAHAHAHAHHAHA",
        'ligandName': 'LALSDLASDLASLDASLDA IAM CRAZYYYY',
        'pdb': open('test.pdb').read(),
        'sdf': open('test.sdf').read(),
        'affinity': 10
    },
        {
            'experimentName': 'Experiment 2',
            'proteinName': "AHAHAHAHAHHAHA",
            'ligandName': 'LALSDLASDLASLDASLDA IAM CRAZYYYY',
            'pdb': open('test.pdb').read(),
            'sdf': open('test.sdf').read(),
            'affinity': 10
        }
    ]
    return experiments