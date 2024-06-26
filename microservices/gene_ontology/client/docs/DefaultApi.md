# gene_ontology_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**is_job_running_job_job_id_is_running_get**](DefaultApi.md#is_job_running_job_job_id_is_running_get) | **GET** /job/{job_id}/is-running | Is Job Running
[**run_go_prediction_run_post**](DefaultApi.md#run_go_prediction_run_post) | **POST** /run | Run Go Prediction


# **is_job_running_job_job_id_is_running_get**
> IsJobRunningResponse is_job_running_job_job_id_is_running_get(job_id)

Is Job Running

### Example


```python
import time
import os
import gene_ontology_microservice
from gene_ontology_microservice.models.is_job_running_response import IsJobRunningResponse
from gene_ontology_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = gene_ontology_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with gene_ontology_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = gene_ontology_microservice.DefaultApi(api_client)
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

# **run_go_prediction_run_post**
> RunGeneOntologyPredictionResponse run_go_prediction_run_post(run_gene_ontology_prediction_request)

Run Go Prediction

### Example


```python
import time
import os
import gene_ontology_microservice
from gene_ontology_microservice.models.run_gene_ontology_prediction_request import RunGeneOntologyPredictionRequest
from gene_ontology_microservice.models.run_gene_ontology_prediction_response import RunGeneOntologyPredictionResponse
from gene_ontology_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = gene_ontology_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with gene_ontology_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = gene_ontology_microservice.DefaultApi(api_client)
    run_gene_ontology_prediction_request = gene_ontology_microservice.RunGeneOntologyPredictionRequest() # RunGeneOntologyPredictionRequest | 

    try:
        # Run Go Prediction
        api_response = api_instance.run_go_prediction_run_post(run_gene_ontology_prediction_request)
        print("The response of DefaultApi->run_go_prediction_run_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->run_go_prediction_run_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **run_gene_ontology_prediction_request** | [**RunGeneOntologyPredictionRequest**](RunGeneOntologyPredictionRequest.md)|  | 

### Return type

[**RunGeneOntologyPredictionResponse**](RunGeneOntologyPredictionResponse.md)

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

