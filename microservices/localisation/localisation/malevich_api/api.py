import pandas as pd

from malevich.square import processor, DF

from localisation.loggers import logger
from localisation.malevich_api.api_models import RunLocalisationPredictionResponse, RunLocalisationPredictionRequest
from localisation.services import run_localisation


@processor()
async def run_localisation_prediction(request: DF[RunLocalisationPredictionRequest]):
    model = request.iat[0].map_to()

    logger.run_localisation_prediction_request(model)
    result = run_localisation(model)
    logger.run_localisation_prediction_response(model)
    return pd.DataFrame([RunLocalisationPredictionResponse.map_from(result)])
