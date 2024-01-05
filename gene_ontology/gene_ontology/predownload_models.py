from esm import pretrained  # type: ignore
from gene_ontology.settings import Settings

settings = Settings()

pretrained.load_model_and_alphabet(settings.model_weights_name)
