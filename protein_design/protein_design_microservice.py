import argparse

from flask_cors import CORS
from flask import request, jsonify, Flask

app = Flask(__name__)
CORS(app)


@app.route('/conformations', methods=['POST'])
def conformations():
    data = request.get_json()
    return jsonify({'protein_design_result': simulation_result})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', default=5785)
    parser.add_argument('--host', default='127.0.0.1')
    args = vars(parser.parse_args())
    host = args['host']
    port = args['port']

    print('--Starting protein design microservice')
    app.run(debug=True, host=host, port=port)