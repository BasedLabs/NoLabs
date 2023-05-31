from flask import Flask, request, jsonify
import os

from .config_reader import get_models_from_config
from src.ai.model_factory import create_model

app = Flask(__name__)

models_metadata = get_models_from_config()

@app.route("/", methods=["POST"])
def index():
    if request.method == "POST":
        response = {}
        sequence = request.json["sequence"]

        models_metadata = get_models_from_config()

        for model_metadata in models_metadata:
            model = create_model(model_metadata)
            response[model.model_name] = model.predict(sequence)

        return jsonify(response)

if __name__ == "__main__":
    app.run(debug=True)
