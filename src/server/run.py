import argparse

from flask import Flask, render_template, request
import src.server.services.inference_service as inference_service
from src.server.services.oboreader import read_obo
from src.server.services.fasta_reader import get_sequences

app = Flask(__name__)
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.config['FLASK_DEBUG'] = True
parser = argparse.ArgumentParser()
parser.add_argument('--port', default=5000)
parser.add_argument('--host', default='127.0.0.1')
parser.add_argument('--test', default=False)
is_test = False


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/amino-acid-page')
def amino_acid_page():
    return render_template('amino-acid-page.html')


@app.route('/drug-target-discovery')
def drug_target_discovery():
    return render_template('drug-target-discovery.html')


@app.route('/inference-amino-acid-page', methods=['POST'])
def inference_amino_acid():
    amino_acid_input_sequence = request.form['inputSequence']
    amino_acid_input_sequence_files = request.files.getlist('inputSequenceFile')
    if not amino_acid_input_sequence:
        amino_acid_input_sequence = \
            [seq for seq in get_sequences(amino_acid_input_sequence_files)][0]
    pipeline = inference_service.create_pipeline(use_gpu=True, is_test=is_test)
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


@app.route('/inference-drug-target-discovery', methods=['POST'])
def inference_drug_target_discovery():
    ligand_files = request.files.getlist('smilesFileInput')
    protein_files = request.files.getlist('proteinFileInput')

    pipeline = inference_service.create_pipeline(is_test=is_test)

    inference_service.save_uploaded_files(pipeline, ligand_files)
    inference_service.save_uploaded_files(pipeline, protein_files)


    inference_service.generate_dti_results(pipeline, ligand_files, protein_files)
    pdb_content, protein_name, ligands_sdf_contents, ligand_names, affinity_list \
         = inference_service.get_dti_results(pipeline, ligand_files)

    if not ligand_files:
        return 'You must provide either smiles text or smiles files', 400

    if not protein_files:
        return 'You must provide a protein .pdb file', 400

    res = [{'proteinName': protein_name,
            'ligandName': ligand_name,
            'pdb': pdb_content,
            'sdf': ligand_content,
            'affinity': affinity} for ligand_name, ligand_content, affinity \
            in zip(ligand_names, ligands_sdf_contents, affinity_list)]

    return {'drugTarget': res}


if __name__ == '__main__':
    args = vars(parser.parse_args())
    is_test = args['test']
    app.run(host=args['host'], port=args['port'])
