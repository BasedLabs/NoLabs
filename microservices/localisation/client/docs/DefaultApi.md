# localisation_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**is_job_running_job_job_id_is_running_get**](DefaultApi.md#is_job_running_job_job_id_is_running_get) | **GET** /job/{job_id}/is-running | Is Job Running
[**predict_run_post**](DefaultApi.md#predict_run_post) | **POST** /run | Predict


# **is_job_running_job_job_id_is_running_get**
> IsJobRunningResponse is_job_running_job_job_id_is_running_get(job_id)

Is Job Running

### Example


```python
import time
import os
import localisation_microservice
from localisation_microservice.models.is_job_running_response import IsJobRunningResponse
from localisation_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = localisation_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with localisation_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = localisation_microservice.DefaultApi(api_client)
    job_id = 'job_id_example' # str | 

    try:
        # Is Job Running
        api_response = api_instance.is_job_running_job_job_id_is_running_get(job_id)
        print("The response of DefaultApi->is_job_running_job_job_id_is_running_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->is_job_running_job_job_id_is_running_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **job_id** | **str**|  | 

### Return type

[**IsJobRunningResponse**](IsJobRunningResponse.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details

| Status code | Description | Response headers |
|-------------|-------------|------------------|
**200** | Successful Response |  -  |
**422** | Validation Error |  -  |

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

# **predict_run_post**
> RunLocalisationPredictionResponse predict_run_post(run_localisation_prediction_request)

Predict

### Example


```python
import time
import os
import localisation_microservice
from localisation_microservice.models.run_localisation_prediction_request import RunLocalisationPredictionRequest
from localisation_microservice.models.run_localisation_prediction_response import RunLocalisationPredictionResponse
from localisation_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = localisation_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with localisation_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = localisation_microservice.DefaultApi(api_client)
    run_localisation_prediction_request = localisation_microservice.RunLocalisationPredictionRequest() # RunLocalisationPredictionRequest | 

    try:
        # Predict
        api_response = api_instance.predict_run_post(run_localisation_prediction_request)
        print("The response of DefaultApi->predict_run_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->predict_run_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **run_localisation_prediction_request** | [**RunLocalisationPredictionRequest**](RunLocalisationPredictionRequest.md)|  | 

### Return type

[**RunLocalisationPredictionResponse**](RunLocalisationPredictionResponse.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: application/json
 - **Accept**: application/json

### HTTP response details

| Status code | Description | Response headers |
|-------------|-------------|------------------|
**200** | Successful Response |  -  |
**422** | Validation Error |  -  |

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

