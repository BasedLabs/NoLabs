import json
import argparse

import torch
from flask import Flask, render_template, request
import src.server.services.inference_service as inference_service
from src.server.services.fasta_reader import get_sequences

app = Flask(__name__)
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.config['FLASK_DEBUG'] = True
parser = argparse.ArgumentParser()
parser.add_argument('--port', default=5000)
parser.add_argument('--host', default='127.0.0.1')


@app.route('/')
def hello():
    return render_template('index.html')


@app.route('/inference', methods=['POST'])
def inference():
    amino_acid_input_sequence = request.form['inputSequence']
    amino_acid_input_sequence_files = request.files.getlist('inputSequenceFile')
    if not amino_acid_input_sequence and not amino_acid_input_sequence_files:
        return render_template('error.html', error='You must enter sequence or select a FASTA file')
    if not amino_acid_input_sequence:
        amino_acid_input_sequence = \
        [seq for seq in get_sequences(amino_acid_input_sequence_files)][0]
    gpu = request.form.get('gpu', False)
    if gpu and not torch.cuda.is_available():
        return render_template('error.html', error="You don't have a CUDA installed or NVIDIA videocard")
    # aminoAcidInputSequenceFile = request.files['inputSequenceFile'] # will consider this later
    pipeline = inference_service.create_pipeline(use_gpu=gpu)
    localisation_result = inference_service.get_localisation_output(pipeline=pipeline,
                                                                    amino_acid_sequence=amino_acid_input_sequence)

    folding_result = inference_service.get_folding_output(pipeline=pipeline,
                                                          amino_acid_sequence=amino_acid_input_sequence)

    esm_protein_localization = {key: value for key, value in localisation_result}

    inference_data = json.dumps({
        'sequence': amino_acid_input_sequence,
        'localisation': {
            'mithochondria': esm_protein_localization["Mitochondrial Proteins"],
            'nucleus': esm_protein_localization["Nuclear Proteins"],
            'cytoplasm': esm_protein_localization["Cytosolic Proteins"],
            'other': esm_protein_localization["Other Proteins"],
            'extracellular': esm_protein_localization["Extracellular/Secreted Proteins"]
        },
        'folding': folding_result
    })
    return render_template('index.html', inference=inference_data)


if __name__ == '__main__':
    args = vars(parser.parse_args())
    app.run(host=args['host'], port=args['port'])
