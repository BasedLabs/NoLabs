from esmfold_light.model import APIFolding
from esmfold_light.api_models import RunEsmFoldPredictionRequest, RunEsmFoldPredictionResponse
from esmfold_light.loggers import Log

__all__ = ['run_facebook_api_folding']


def run_facebook_api_folding(parameters: RunEsmFoldPredictionRequest) -> RunEsmFoldPredictionResponse:
    try:
        model = APIFolding()
        pdb_content = model.predict(
            sequence=parameters.protein_sequence
        )
        return RunEsmFoldPredictionResponse(pdb_content=pdb_content)
    except Exception as e:
        Log.exception()
        return RunEsmFoldPredictionResponse(pdb_content=None,
                                            errors=['Unable to run esmfold. Internal error', str(e)])
