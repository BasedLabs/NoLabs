from flask import Request

from src.server.services.conformations.conformations_pipeline import pipeline
from src.server.services.experiments_structure_loader import ConformationsExperimentsLoader
from src.server.api_handlers.api_handler import ApiHandler


class ConformationsApiHandler(ApiHandler):
    def __init__(self):
        self.experiments_loader = ConformationsExperimentsLoader()

    def get_experiments(self):
        return self.experiments_loader.load_experiments()

    def get_experiment(self, request):
        experiment_id = request.args.get('id')
        experiment_name = request.args.get('name')

        if not self.experiments_loader.experiment_exists(experiment_id):
            return {'id': experiment_id, 'name': experiment_name, 'data': {}}

        experiment_data = self.experiments_loader.load_experiment(experiment_id)
        return {'id': experiment_id, 'name': experiment_name, 'data': experiment_data}

    def change_experiment_name(self, request: Request):
        j = request.get_json(force=True)
        experiment_id = j['id']
        experiment_name = j['name']

        self.experiments_loader.rename_experiment(experiment_id, experiment_name)

        return {'status': 200}

    def delete_experiment(self, request: Request):
        j = request.get_json(force=True)
        self.experiments_loader.delete_experiment(j['id'])

        return {'status': 200}

    def inference(self, request: Request) -> dict:
        protein_files = request.files.getlist('proteinFileInput')
        experiment_name = request.form['experimentName']
        experiment_id = request.form['experimentId']
        total_frames = int(request.form['totalFramesInput'])
        system_temp = float(request.form['tempInput'])
        take_frame_every = int(request.form['takeFrameEveryInput'])
        step_size = float(request.form['stepSizeInput'])
        replace_nonstandard_residues = 'replaceNonstandardResiduesInput' in request.form
        integrator = request.form['integratorSelect']
        friction_coeff = float(request.form['frictionCoeffInput'])
        add_missing_hydrogens = 'addMissingHydrogensCheckbox' in request.form
        add_missing_atoms = 'addMissingAtomsCheckbox' in request.form
        ignore_missing_atoms = 'ignoreMissingAtomsCheckbox' in request.form

        simulation_result = pipeline(pdb_content=protein_files[0],
                                     total_frames=total_frames,
                                     take_frame_every=take_frame_every,
                                     integrator=integrator,
                                     friction_coeff=friction_coeff,
                                     step_size=step_size,
                                     temperature_kelvin=system_temp,
                                     replace_nonstandard_residues=replace_nonstandard_residues,
                                     add_missing_atoms=add_missing_atoms,
                                     add_missing_hydrogens=add_missing_hydrogens,
                                     ignore_missing=ignore_missing_atoms)
        if simulation_result:
            self.experiments_loader.store_experiment(experiment_id, simulation_result)
            self.experiments_loader.save_experiment_metadata(experiment_id, experiment_name)

        return {
            'id': experiment_id,
            'name': experiment_name,
            'data': {'pdb': simulation_result}
        }
