# import ssl
# ssl._create_default_https_context = ssl._create_unverified_context
import argparse

import torch

from nolabs.server import settings

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', default=5000)
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--test', action='store_true')
    parser.add_argument('--demo', action='store_true')
    args = vars(parser.parse_args())
    is_test = args['test']
    settings.is_test = is_test
    is_demo = args['demo']
    settings.is_demo = is_demo
    settings.use_gpu = torch.cuda.is_available()
    settings.host = args['host']
    settings.port = args['port']

    from nolabs.server.initializers.initialize_app import init

    app, socketio = init()

    from nolabs.server import factories
    import nolabs.server.amino_acid_blueprint as amino_acid_blueprint
    import nolabs.server.drug_target_blueprint as drug_target_blueprint
    import nolabs.server.conformations_blueprint as conformations_blueprint
    import nolabs.server.protein_design_blueprint as protein_design_blueprint
    import nolabs.server.logging_blueprint as logging_blueprint

    amino_acid_api_handler, drug_target_api_handler, conformations_api_handler, protein_design_api_handler \
        = factories.api_handlers_factory(is_test=is_test, is_demo=is_demo)

    app.register_blueprint(drug_target_blueprint.resolve_api_endpoints(drug_target_api_handler),
                           url_prefix='/api/drug-target')
    app.register_blueprint(amino_acid_blueprint.resolve_api_endpoints(amino_acid_api_handler),
                           url_prefix='/api/amino-acid')
    app.register_blueprint(conformations_blueprint.resolve_api_endpoints(conformations_api_handler),
                           url_prefix='/api/conformations')
    app.register_blueprint(protein_design_blueprint.resolve_api_endpoints(protein_design_api_handler),
                           url_prefix='/api/protein-design')
    app.register_blueprint(logging_blueprint.resolve_api_endpoints(), url_prefix='/logging')
    app.config['TEMPLATES_AUTO_RELOAD'] = True
    app.config['FLASK_DEBUG'] = True
    print('-- Starting flask server')
    socketio.run(app, host=settings.host, port=settings.port, allow_unsafe_werkzeug=True)
