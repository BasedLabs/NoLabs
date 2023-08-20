import json
import argparse

import torch
from flask import Flask, render_template, request
import src.server.services.inference_service as inference_service
from server.services.oboreader import read_obo
from server.services.obo_loader_service import load_obo
from src.server.services.fasta_reader import get_sequences

app = Flask(__name__)
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.config['FLASK_DEBUG'] = True
parser = argparse.ArgumentParser()
parser.add_argument('--port', default=5005)
parser.add_argument('--host', default='127.0.0.1')
parser.add_argument('--test', default=False)
is_test = False


@app.route('/')
def hello():
    return render_template('index.html')


@app.route('/obo-graph-data')
def cytoscape_data():
    obo_graph, obo_graph_light = read_obo()
    return {
        'obo_graph': obo_graph,
        'obo_graph_light': obo_graph_light
    }

@app.route('/cytoscape')
def cytoscape():
    return render_template('cytoscape.html')

@app.route('/d3')
def d3():
    return render_template('d3.html')


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
    pipeline = inference_service.create_pipeline(use_gpu=gpu, is_test=is_test)
    localisation_result = inference_service.get_localisation_output(pipeline=pipeline,
                                                                    amino_acid_sequence=amino_acid_input_sequence)

    folding_result = inference_service.get_folding_output(pipeline=pipeline,
                                                          amino_acid_sequence=amino_acid_input_sequence)

    esm_protein_localization = {key: value for key, value in localisation_result}

    obo_ids = {
        "GO:0009274",
        "GO:0000724",
        "GO:0009314",
        "GO:0009432",
        "GO:0005524",
        "GO:0006281",
        "GO:0006310",
        "GO:0005524",
        "GO:0006281",
        "GO:0006974",
        "GO:0000166"
    }

    inference_data = json.dumps({
        'sequence': amino_acid_input_sequence,
        'localisation': {
            'mithochondria': esm_protein_localization["Mitochondrial Proteins"],
            'nucleus': esm_protein_localization["Nuclear Proteins"],
            'cytoplasm': esm_protein_localization["Cytosolic Proteins"],
            'other': esm_protein_localization["Other Proteins"],
            'extracellular': esm_protein_localization["Extracellular/Secreted Proteins"]
        },
        'folding': folding_result,
        'gene_ontology': load_obo(ids=obo_ids)
    })
    return render_template('index.html', inference=inference_data)


if __name__ == '__main__':
    args = vars(parser.parse_args())
    is_test = args['test']
    app.run(host=args['host'], port=args['port'])
