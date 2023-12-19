import argparse

from flask_cors import CORS
from flask import request, jsonify, Flask
import conformations_pipeline as pipeline

app = Flask(__name__)
CORS(app)


@app.route('/conformations', methods=['POST'])
def conformations():
    data = request.get_json()
    pdb_content = data.get('pdb_content')
    total_frames = data.get('total_frames', 10000)
    take_frame_every = data.get('take_frame_every', 1000)
    integrator = data.get('integrator', 'LangevinIntegator')
    friction_coeff = data.get('friction_coeff', 1.0)
    step_size = data.get('step_size', 0.002)
    temperature_kelvin = data.get('temperature_kelvin', 273.15)
    replace_nonstandard_residues = data.get('replace_nonstandard_residues', True)
    add_missing_atoms = data.get('add_missing_atoms', True)
    add_missing_hydrogens = data.get('add_missing_hydrogens', True)
    ignore_missing = data.get('ignore_missing', False)
    simulation_result = pipeline.pipeline(pdb_content,
                                          total_frames,
                                          take_frame_every,
                                          integrator,
                                          friction_coeff,
                                          step_size,
                                          temperature_kelvin,
                                          replace_nonstandard_residues,
                                          add_missing_atoms,
                                          add_missing_hydrogens,
                                          ignore_missing)
    return jsonify({'simulation_result': simulation_result})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', default=5785)
    parser.add_argument('--host', default='127.0.0.1')
    args = vars(parser.parse_args())
    host = args['host']
    port = args['port']

    print('--Starting conformations microservice')
    app.run(debug=True, host=host, port=port)
