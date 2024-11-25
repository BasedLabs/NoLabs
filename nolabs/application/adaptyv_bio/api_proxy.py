from typing import Literal, Dict, Any, List, Optional

import requests
from pydantic import BaseModel, Field
from requests.adapters import HTTPAdapter, Retry

from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.infrastructure.log import logger
from nolabs.infrastructure.settings import settings, Environment


class Target(BaseModel):
    id: str
    name: str
    description: str
    swissprot_id: str


class GetExperimentEstimatesResponse(BaseModel):
    total_price: int = Field(..., alias='totalPrice')
    turnaround_time: int = Field(..., alias='turnaroundTime')


class GetExperimentsTypesResponse(BaseModel):
    id: str
    full_name: str


class AdaptyvBioApi:
    def _make_request(self,
                      method: Literal['GET', 'POST'],
                      endpoint: str,
                      json: Optional[Dict[str, Any]] = None,
                      params: Optional[Dict[str, Any]] = None,
                      expect_response=True):
        url = f"{settings.adaptyv_bio_api_base}/{endpoint}"
        headers = {"Authorization": f"Bearer {settings.adaptyv_bio_api_token}"}

        try:
            s = requests.Session()
            retries = Retry(total=5, backoff_factor=1, status_forcelist=[502, 503, 504, 429])
            s.mount('https://', HTTPAdapter(max_retries=retries))

            response = s.request(
                method=method,
                url=url,
                headers=headers,
                json=json,
                params=params,
                timeout=10.0
            )
            response.raise_for_status()
            if expect_response:
                return response.json()
        except requests.exceptions.HTTPError as http_err:
            if http_err.response.status_code in {401, 403}:
                raise NoLabsException(ErrorCodes.adaptyv_bio_api_unauthorized)
            if http_err.response.status_code == 429:
                raise NoLabsException(ErrorCodes.adaptyv_bio_api_too_many_requests)
            logger.exception("Adaptyv bio exception occurred")
            raise NoLabsException(ErrorCodes.adaptyv_bio_api_error)
        except requests.exceptions.RequestException as e:
            logger.exception("Adaptyv bio exception occurred")
            raise NoLabsException(ErrorCodes.adaptyv_bio_api_error)

    def _get_experiments_types(self) -> List[GetExperimentsTypesResponse]:
        response = self._make_request('GET', 'get_experiment_types')
        return [GetExperimentsTypesResponse(**j) for j in response]

    def list_targets(self, query: str) -> List[Target]:
        response = self._make_request('GET', 'list_targets', params={'searchQuery': query})
        return [Target(**p) for p in response]


class AdaptyvBioProteinBindingScreeningApi(AdaptyvBioApi):
    experiment_type_id: str = '45b59ab3-2402-4001-a401-b3e1203c5ae5'

    def get_experiment_estimates(self, n_designs: int, avg_length: int, n_replicates: int):
        response = self._make_request('POST', 'get_experiment_estimates',
                                      json={
                                          'experiment_type_id': self.experiment_type_id,
                                          'params': {
                                              'n_designs': n_designs,
                                              'avg_length': avg_length,
                                              'n_replicates': n_replicates,
                                          }
                                      })
        return GetExperimentEstimatesResponse(**response)

    def submit_experiment(self,
                          sequences: List[str],
                          target_id: str,
                          email: str,
                          session_url: str,
                          n_replicates: int,
                          cart_total: int,
                          avg_length: int,
                          n_designs: int
                          ):
        json_data = {
            'sequences': sequences,
            'target_id': target_id,
            'email': email,
            'experiment_type_id': self.experiment_type_id,
            'session_url': session_url,
            'n_replicates': n_replicates,
            'cart_total': cart_total,
            'data': {
                'avg_length': avg_length,
                'n_designs': n_designs
            }
        }

        if settings.environment != Environment.production:
            json_data['data']['env'] = 'test'

        self._make_request(method='POST',
                           endpoint='submit_experiment',
                           json={
                               'sequences': sequences,
                               'target_id': target_id,
                               'email': email,
                               'experiment_type_id': self.experiment_type_id,
                               'session_url': session_url,
                               'n_replicates': n_replicates,
                               'cart_total': cart_total,
                               'data': {
                                   'avg_length': avg_length,
                                   'n_designs': n_designs
                               }
                           },
                           expect_response=False)



class AdaptyvBioProteinAffinityCharacterizationApi(AdaptyvBioApi):
    experiment_type_id: str = '2c17531c-5a8b-45a7-9b3f-5f82e34c050d'

    def get_experiment_estimates(self, n_designs: int, avg_length: int, n_replicates: int):
        response = self._make_request('POST', 'get_experiment_estimates',
                                      json={
                                          'experiment_type_id': self.experiment_type_id,
                                          'params': {
                                              'n_designs': n_designs,
                                              'avg_length': avg_length,
                                              'n_replicates': n_replicates,
                                          }
                                      })
        return GetExperimentEstimatesResponse(**response)

    def submit_experiment(self,
                          sequences: List[str],
                          target_id: str,
                          email: str,
                          session_url: str,
                          n_replicates: int,
                          cart_total: int,
                          avg_length: int,
                          n_designs: int
                          ):
        json_data = {
            'sequences': sequences,
            'target_id': target_id,
            'email': email,
            'experiment_type_id': self.experiment_type_id,
            'session_url': session_url,
            'n_replicates': n_replicates,
            'cart_total': cart_total,
            'data': {
                'avg_length': avg_length,
                'n_designs': n_designs
            }
        }

        if settings.environment != Environment.production:
            json_data['data']['env'] = 'test'

        self._make_request(method='POST',
                           endpoint='submit_experiment',
                           json={
                               'sequences': sequences,
                               'target_id': target_id,
                               'email': email,
                               'experiment_type_id': self.experiment_type_id,
                               'session_url': session_url,
                               'n_replicates': n_replicates,
                               'cart_total': cart_total,
                               'data': {
                                   'avg_length': avg_length,
                                   'n_designs': n_designs
                               }
                           },
                           expect_response=False)
