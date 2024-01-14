# localisation_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**run_localisation_prediction_run_localisation_prediction_post**](DefaultApi.md#run_localisation_prediction_run_localisation_prediction_post) | **POST** /run-localisation-prediction | Run Localisation Prediction


# **run_localisation_prediction_run_localisation_prediction_post**
> RunLocalisationPredictionResponse run_localisation_prediction_run_localisation_prediction_post(run_localisation_prediction_request)

Run Localisation Prediction

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
        # Run Localisation Prediction
        api_response = api_instance.run_localisation_prediction_run_localisation_prediction_post(run_localisation_prediction_request)
        print("The response of DefaultApi->run_localisation_prediction_run_localisation_prediction_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->run_localisation_prediction_run_localisation_prediction_post: %s\n" % e)
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

