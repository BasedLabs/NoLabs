from esmfold.model import Folding, APIFolding
from esmfold.api_models import RunEsmFoldPredictionRequest, RunEsmFoldPredictionResponce
from esmfold.loggers import Log

__all__ = ['run_folding', 'run_facebook_api_folding']

def run_folding(parameters: RunEsmFoldPredictionRequest) -> RunEsmFoldPredictionResponce:
    try:
        model = Folding()
        model.load_model()
        pdb_content = model.predict(
            sequence=parameters.protein_sequence
        )
        return RunEsmFoldPredictionResponce(pdb_content=pdb_content)
    except Exception as e:
        Log.exception()
        return RunEsmFoldPredictionResponce(pdb_content=None,
                                         errors=['Unable to run esmfold. Internal error', str(e)])


def run_facebook_api_folding(parameters: RunEsmFoldPredictionRequest) -> RunEsmFoldPredictionResponce:
    try:
        model = APIFolding()
        pdb_content = model.predict(
            sequence=parameters.protein_sequence
        )
        return RunEsmFoldPredictionResponce(pdb_content=pdb_content)
    except Exception as e:
        Log.exception()
        return RunEsmFoldPredictionResponce(pdb_content=None,
                                            errors=['Unable to run esmfold. Internal error', str(e)])
