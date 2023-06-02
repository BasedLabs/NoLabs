from src.ai import model_factory
from src.ai.model import BaseModel


def get_localisation_output(amino_acid_sequence: str) -> str:
    assert amino_acid_sequence
    localisation_model = model_factory.create_model({'name': 'cell localisation response',
                                                    'type': 'classification'})
    return localisation_model.predict(amino_acid_sequence)