from esmfold.model import Folding
from esmfold.api_models import RunEsmFoldPredictionRequest, RunEsmFoldPredictionResponse
from esmfold.loggers import Log
import torch

__all__ = ['run_folding']

Log.cuda_is_avaialable(torch.cuda.is_available())

def run_folding(parameters: RunEsmFoldPredictionRequest) -> RunEsmFoldPredictionResponse:
    try:
        model = Folding()
        model.load_model()
        pdb_content = model.predict(
            sequence=parameters.protein_sequence
        )
        return RunEsmFoldPredictionResponse(pdb_content=pdb_content, errors=[])
    except Exception as e:
        Log.exception()
        return RunEsmFoldPredictionResponse(pdb_content=None,
                                         errors=['Unable to run esmfold. Internal error', str(e)])
