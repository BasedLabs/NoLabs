import argparse

from flask_cors import CORS
from flask import request, jsonify, Flask
from protein_design_pipeline import pipeline

app = Flask(__name__)
CORS(app)


@app.route('/protein-design', methods=['POST'])
def protein_design():
    data = request.get_json()

    pdb_content = data.get('pdb_content', None)
    contig = data.get('contig', '50')
    symmetry = data.get('symmetry', None)
    timesteps = data.get('timesteps', None)
    hotspots = data.get('hotspots', None)

    inference_result = pipeline(pdb_content, contig, symmetry, timesteps, hotspots)

    return jsonify({'result': inference_result})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', default=5786)
    parser.add_argument('--host', default='127.0.0.1')
    args = vars(parser.parse_args())
    host = args['host']
    port = args['port']

    print('--Starting protein design microservice')
    app.run(debug=True, host=host, port=port)