from flask import Blueprint

from flask import Flask, render_template, request
import src.server.services.inference_service as inference_service
from src.server import settings
from src.server.services.oboreader import read_obo
from src.server.services.fasta_reader import get_sequences

amino_acid_bp = Blueprint('amino-acid', __name__,
                          template_folder='templates')


@amino_acid_bp.route('/')
def amino_acid_page():
    return render_template('amino-acid-page.html')


@amino_acid_bp.route('/amino-acid-inference', methods=['POST'])
def inference_amino_acid():
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
