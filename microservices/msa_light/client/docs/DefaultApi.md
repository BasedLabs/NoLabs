# msa_light_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**get_running_jobs_jobs_running_get**](DefaultApi.md#get_running_jobs_jobs_running_get) | **GET** /jobs/running | Get Running Jobs
[**is_job_running_job_job_id_is_running_get**](DefaultApi.md#is_job_running_job_job_id_is_running_get) | **GET** /job/{job_id}/is-running | Is Job Running
[**predict_msa_predict_msa_post**](DefaultApi.md#predict_msa_predict_msa_post) | **POST** /predict-msa | Predict Msa


# **get_running_jobs_jobs_running_get**
> object get_running_jobs_jobs_running_get()

Get Running Jobs

### Example


```python
import time
import os
import msa_light_microservice
from msa_light_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = msa_light_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with msa_light_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = msa_light_microservice.DefaultApi(api_client)

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
> object is_job_running_job_job_id_is_running_get(job_id)

Is Job Running

### Example


```python
import time
import os
import msa_light_microservice
from msa_light_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = msa_light_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with msa_light_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = msa_light_microservice.DefaultApi(api_client)
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
**422** | Validation Error |  -  |

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

# **predict_msa_predict_msa_post**
> RunMsaPredictionResponse predict_msa_predict_msa_post(run_msa_prediction_request)

Predict Msa

### Example


```python
import time
import os
import msa_light_microservice
from msa_light_microservice.models.run_msa_prediction_request import RunMsaPredictionRequest
from msa_light_microservice.models.run_msa_prediction_response import RunMsaPredictionResponse
from msa_light_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = msa_light_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with msa_light_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = msa_light_microservice.DefaultApi(api_client)
    run_msa_prediction_request = msa_light_microservice.RunMsaPredictionRequest() # RunMsaPredictionRequest | 

    try:
        # Predict Msa
        api_response = api_instance.predict_msa_predict_msa_post(run_msa_prediction_request)
        print("The response of DefaultApi->predict_msa_predict_msa_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->predict_msa_predict_msa_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **run_msa_prediction_request** | [**RunMsaPredictionRequest**](RunMsaPredictionRequest.md)|  | 

### Return type

[**RunMsaPredictionResponse**](RunMsaPredictionResponse.md)

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

