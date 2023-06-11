import json
import argparse
import torch
from flask import Flask, render_template, request, jsonify
import src.server.services.inference_service as inference_service

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
    amino_acid_input_sequence_file = request.files.get('inputSequenceFile', None)
    if not amino_acid_input_sequence and not amino_acid_input_sequence_file:
        raise Exception('Input amino acid sequence')
    if not amino_acid_input_sequence:
        amino_acid_input_sequence = amino_acid_input_sequence_file.read().replace('\n', '')
    gpu = request.form.get('gpu', False)
    if gpu and not torch.cuda.is_available():
        return {'result': 'You do not have a gpu installed or CUDA package'}
    # aminoAcidInputSequenceFile = request.files['inputSequenceFile'] # will consider this later
    pipeline = inference_service.create_pipeline()
    localisation_result = inference_service.get_localisation_output(pipeline=pipeline,
                                                                  amino_acid_sequence=amino_acid_input_sequence or amino_acid_input_sequence_file)
    
    print(localisation_result)
    
    folding_result = inference_service.get_folding_output(pipeline=pipeline,
                                                                amino_acid_sequence=amino_acid_input_sequence or amino_acid_input_sequence_file)

    esm_protein_localization = {key: value for key, value in localisation_result}

    inference_data = json.dumps({
        'sequence': amino_acid_input_sequence or amino_acid_input_sequence_file,
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