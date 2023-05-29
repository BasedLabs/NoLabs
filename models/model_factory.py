from .model import LocalisationModel

def create_model(model_metadata):
    model_name = model_metadata.get('name')
    model_type = model_metadata.get('type')

    if model_type == "localisation":
        model_labels = model_metadata.get('labels', [])

        model = LocalisationModel(model_name)
        model.load_model()
        model.set_labels(model_labels)

        return model