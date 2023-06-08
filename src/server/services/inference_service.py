from src.ai import model_factory
from src.ai.model import BaseModel
from src.ai.pipeline import Pipeline


def create_pipeline():
    pipeline = Pipeline()
    models_metadata = get_models_from_config()

    for model_metadata in models_metadata:
        model = create_model(model_metadata)
        pipeline.add_model(model)

    return pipeline


def get_localisation_output(pipeline, amino_acid_sequence: str) -> str:
    assert amino_acid_sequence
    model = pipeline.get_model_by_task("localisation")
    return model.predict(amino_acid_sequence)

def read_json_file(file_path):
    with open(file_path, 'r') as json_file:
        data = json.load(json_file)
    return data

def read_config():
    dirname=os.path.dirname
    root_directory = dirname(dirname(os.path.abspath(__file__)))
    file_path = os.path.join(root_directory, 'config.json')
    data = read_json_file(file_path)
    return data

def get_models_from_config():
    data = read_config()
    models = data.get('models', [])
    return models