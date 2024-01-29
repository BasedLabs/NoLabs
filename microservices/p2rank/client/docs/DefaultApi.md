# p2rank_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**get_running_jobs_jobs_running_get**](DefaultApi.md#get_running_jobs_jobs_running_get) | **GET** /jobs/running | Get Running Jobs
[**is_job_running_job_job_id_is_running_get**](DefaultApi.md#is_job_running_job_job_id_is_running_get) | **GET** /job/{job_id}/is-running | Is Job Running
[**predict_run_p2rank_post**](DefaultApi.md#predict_run_p2rank_post) | **POST** /run-p2rank | Predict


# **get_running_jobs_jobs_running_get**
> object get_running_jobs_jobs_running_get()

Get Running Jobs

### Example


```python
import time
import os
import p2rank_microservice
from p2rank_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = p2rank_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with p2rank_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = p2rank_microservice.DefaultApi(api_client)

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
import p2rank_microservice
from p2rank_microservice.models.is_job_running_response import IsJobRunningResponse
from p2rank_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = p2rank_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with p2rank_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = p2rank_microservice.DefaultApi(api_client)
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

# **predict_run_p2rank_post**
> RunP2RankPredictionResponse predict_run_p2rank_post(run_p2_rank_prediction_request)

Predict

### Example


```python
import time
import os
import p2rank_microservice
from p2rank_microservice.models.run_p2_rank_prediction_request import RunP2RankPredictionRequest
from p2rank_microservice.models.run_p2_rank_prediction_response import RunP2RankPredictionResponse
from p2rank_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = p2rank_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with p2rank_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = p2rank_microservice.DefaultApi(api_client)
    run_p2_rank_prediction_request = p2rank_microservice.RunP2RankPredictionRequest() # RunP2RankPredictionRequest | 

    try:
        # Predict
        api_response = api_instance.predict_run_p2rank_post(run_p2_rank_prediction_request)
        print("The response of DefaultApi->predict_run_p2rank_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->predict_run_p2rank_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **run_p2_rank_prediction_request** | [**RunP2RankPredictionRequest**](RunP2RankPredictionRequest.md)|  | 

### Return type

[**RunP2RankPredictionResponse**](RunP2RankPredictionResponse.md)

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

