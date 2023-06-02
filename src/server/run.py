import torch
from flask import Flask, render_template, request

app = Flask(__name__)


@app.route('/')
def hello():
    return render_template('index.html')


@app.route('/inference', methods=['POST'])
def inference():
    aminoAcidInputSequence = request.form['inputSequenceFile']
    gpu = request.form['gpu']
    if not torch.cuda.is_available() and gpu:
        return {'result': 'You do not have a gpu installed or CUDA package'}
    # aminoAcidInputSequenceFile = request.files['inputSequenceFile'] # will consider this later

    return {'result': 'ok'}
