# umol_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**predict_run_folding_post**](DefaultApi.md#predict_run_folding_post) | **POST** /run-folding | Predict
[**predict_through_api_run_through_api_folding_post**](DefaultApi.md#predict_through_api_run_through_api_folding_post) | **POST** /run-through-api-folding | Predict Through Api


# **predict_run_folding_post**
> RunEsmFoldPredictionResponse predict_run_folding_post(run_esm_fold_prediction_request)

Predict

### Example


```python
import time
import os
import umol_microservice
from umol_microservice.models.run_esm_fold_prediction_request import RunEsmFoldPredictionRequest
from umol_microservice.models.run_esm_fold_prediction_response import RunEsmFoldPredictionResponse
from umol_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = umol_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with umol_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = umol_microservice.DefaultApi(api_client)
    run_esm_fold_prediction_request = umol_microservice.RunEsmFoldPredictionRequest() # RunEsmFoldPredictionRequest | 

    try:
        # Predict
        api_response = api_instance.predict_run_folding_post(run_esm_fold_prediction_request)
        print("The response of DefaultApi->predict_run_folding_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->predict_run_folding_post: %s\n" % e)
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

# **predict_through_api_run_through_api_folding_post**
> RunEsmFoldPredictionResponse predict_through_api_run_through_api_folding_post(run_esm_fold_prediction_request)

Predict Through Api

### Example


```python
import time
import os
import umol_microservice
from umol_microservice.models.run_esm_fold_prediction_request import RunEsmFoldPredictionRequest
from umol_microservice.models.run_esm_fold_prediction_response import RunEsmFoldPredictionResponse
from umol_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = umol_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with umol_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = umol_microservice.DefaultApi(api_client)
    run_esm_fold_prediction_request = umol_microservice.RunEsmFoldPredictionRequest() # RunEsmFoldPredictionRequest | 

    try:
        # Predict Through Api
        api_response = api_instance.predict_through_api_run_through_api_folding_post(run_esm_fold_prediction_request)
        print("The response of DefaultApi->predict_through_api_run_through_api_folding_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->predict_through_api_run_through_api_folding_post: %s\n" % e)
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

