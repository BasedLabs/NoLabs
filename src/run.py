from flask import Flask, request, jsonify
import json
import plotly
import plotly.express as px
from transformers import AutoTokenizer, AutoModelForSequenceClassification
import torch
import os

app = Flask(__name__)

model_checkpoint_path = 'ritakurban/ESM_protein_localization'
tokenizer = AutoTokenizer.from_pretrained(model_checkpoint_path)
model = AutoModelForSequenceClassification.from_pretrained(model_checkpoint_path)

def predict_localization(sequence):
    inputs = tokenizer(sequence, return_tensors='pt', padding=True, truncation=True)
    with torch.no_grad():
        outputs = model(**inputs)
        probabilities = torch.nn.functional.softmax(outputs.logits, dim=-1).tolist()[0]
    return probabilities

@app.route("/", methods=["POST"])
def index():
    if request.method == "POST":
        print(request.json)
        sequence = request.json["sequence"]
        probabilities = predict_localization(sequence)
        labels = ['Cytosolic Proteins', 'Mitochondrial Proteins', 'Nuclear Proteins', 'Other Proteins', 'Extracellular/Secreted Proteins', ]
        fig = px.bar(x=labels, y=probabilities)
        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        prob_table = list(zip(labels, probabilities))
        return jsonify(prob_table)

if __name__ == "__main__":
    app.run(debug=True)
