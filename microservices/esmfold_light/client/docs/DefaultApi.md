# esmfold_light_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**get_running_jobs_jobs_running_get**](DefaultApi.md#get_running_jobs_jobs_running_get) | **GET** /jobs/running | Get Running Jobs
[**is_job_running_job_job_id_is_running_get**](DefaultApi.md#is_job_running_job_job_id_is_running_get) | **GET** /job/{job_id}/is-running | Is Job Running
[**predict_through_api_run_folding_post**](DefaultApi.md#predict_through_api_run_folding_post) | **POST** /run-folding | Predict Through Api


# **get_running_jobs_jobs_running_get**
> object get_running_jobs_jobs_running_get()

Get Running Jobs

### Example


```python
import time
import os
import esmfold_light_microservice
from esmfold_light_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = esmfold_light_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with esmfold_light_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = esmfold_light_microservice.DefaultApi(api_client)

    try:
        # Get Running Jobs
        api_response = api_instance.get_running_jobs_jobs_running_get()
        print("The response of DefaultApi->get_running_jobs_jobs_running_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->get_running_jobs_jobs_running_get: %s\n" % e)
```



### Parameters

This endpoint does not need any parameter.

### Return type

**object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details

| Status code | Description | Response headers |
|-------------|-------------|------------------|
**200** | Successful Response |  -  |

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

# **is_job_running_job_job_id_is_running_get**
> IsJobRunningResponse is_job_running_job_job_id_is_running_get(job_id)

Is Job Running

### Example


```python
import time
import os
import esmfold_light_microservice
from esmfold_light_microservice.models.is_job_running_response import IsJobRunningResponse
from esmfold_light_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = esmfold_light_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with esmfold_light_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = esmfold_light_microservice.DefaultApi(api_client)
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

# **predict_through_api_run_folding_post**
> RunEsmFoldPredictionResponse predict_through_api_run_folding_post(run_esm_fold_prediction_request)

Predict Through Api

### Example


```python
import time
import os
import esmfold_light_microservice
from esmfold_light_microservice.models.run_esm_fold_prediction_request import RunEsmFoldPredictionRequest
from esmfold_light_microservice.models.run_esm_fold_prediction_response import RunEsmFoldPredictionResponse
from esmfold_light_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = esmfold_light_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with esmfold_light_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = esmfold_light_microservice.DefaultApi(api_client)
    run_esm_fold_prediction_request = esmfold_light_microservice.RunEsmFoldPredictionRequest() # RunEsmFoldPredictionRequest | 

    try:
        # Predict Through Api
        api_response = api_instance.predict_through_api_run_folding_post(run_esm_fold_prediction_request)
        print("The response of DefaultApi->predict_through_api_run_folding_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->predict_through_api_run_folding_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **run_esm_fold_prediction_request** | [**RunEsmFoldPredictionRequest**](RunEsmFoldPredictionRequest.md)|  | 

### Return type

[**RunEsmFoldPredictionResponse**](RunEsmFoldPredictionResponse.md)

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

