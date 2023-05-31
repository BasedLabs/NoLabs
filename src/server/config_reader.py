import json
import os

def read_json_file(file_path):
    with open(file_path, 'r') as json_file:
        data = json.load(json_file)
    return data

def read_config():
    dirname=os.path.dirname
    root_directory = dirname(dirname(dirname(os.path.abspath(__file__))))
    file_path = os.path.join(root_directory, 'config.json')
    data = read_json_file(file_path)
    return data

def get_models_from_config():
    data = read_config()
    models = data.get('models', [])
    return models


