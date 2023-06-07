import json

import torch
from flask import Flask, render_template, request

app = Flask(__name__)
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.config['FLASK_DEBUG'] = True


@app.route('/')
def hello():
    return render_template('index.html')


@app.route('/inference', methods=['POST'])
def inference():
    aminoAcidInputSequence = request.form['inputSequence']
    aminoAcidInputSequenceFile = request.files.get('inputSequenceFile', None)
    gpu = request.form.get('gpu', False)
    if gpu and not torch.cuda.is_available():
        return {'result': 'You do not have a gpu installed or CUDA package'}
    # aminoAcidInputSequenceFile = request.files['inputSequenceFile'] # will consider this later

    inference_data = json.dumps({
            'sequence': aminoAcidInputSequence,
            'localisation': {
                'mithochondria': 0.95,
                'nucleus': 0.01,
                'cytoplasm': 0.01,
            },
            'folding': open('src/server/static/1r6a.pdb', 'r').read()
        })
    return render_template('index.html', inference=inference_data)
