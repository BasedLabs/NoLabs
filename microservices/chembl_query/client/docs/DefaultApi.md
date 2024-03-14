# chembl_query_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**get_running_jobs_jobs_running_get**](DefaultApi.md#get_running_jobs_jobs_running_get) | **GET** /jobs/running | Get Running Jobs
[**is_job_running_job_job_id_is_running_get**](DefaultApi.md#is_job_running_job_job_id_is_running_get) | **GET** /job/{job_id}/is-running | Is Job Running
[**query_query_chembl_by_condition_post**](DefaultApi.md#query_query_chembl_by_condition_post) | **POST** /query-chembl-by-condition | Query
[**query_query_chembl_post**](DefaultApi.md#query_query_chembl_post) | **POST** /query-chembl | Query


# **get_running_jobs_jobs_running_get**
> object get_running_jobs_jobs_running_get()

Get Running Jobs

### Example


```python
import chembl_query_microservice
from chembl_query_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = chembl_query_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with chembl_query_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = chembl_query_microservice.DefaultApi(api_client)

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
import chembl_query_microservice
from chembl_query_microservice.models.is_job_running_response import IsJobRunningResponse
from chembl_query_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = chembl_query_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with chembl_query_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = chembl_query_microservice.DefaultApi(api_client)
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

# **query_query_chembl_by_condition_post**
> DrugIndicationResponse query_query_chembl_by_condition_post(drug_indication_request)

Query

### Example


```python
import chembl_query_microservice
from chembl_query_microservice.models.drug_indication_request import DrugIndicationRequest
from chembl_query_microservice.models.drug_indication_response import DrugIndicationResponse
from chembl_query_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = chembl_query_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with chembl_query_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = chembl_query_microservice.DefaultApi(api_client)
    drug_indication_request = chembl_query_microservice.DrugIndicationRequest() # DrugIndicationRequest | 

    try:
        # Query
        api_response = api_instance.query_query_chembl_by_condition_post(drug_indication_request)
        print("The response of DefaultApi->query_query_chembl_by_condition_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->query_query_chembl_by_condition_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **drug_indication_request** | [**DrugIndicationRequest**](DrugIndicationRequest.md)|  | 

### Return type

[**DrugIndicationResponse**](DrugIndicationResponse.md)

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

# **query_query_chembl_post**
> ChEMBLMoleculeResponse query_query_chembl_post(ch_embl_molecule_request)

Query

### Example


```python
import chembl_query_microservice
from chembl_query_microservice.models.ch_embl_molecule_request import ChEMBLMoleculeRequest
from chembl_query_microservice.models.ch_embl_molecule_response import ChEMBLMoleculeResponse
from chembl_query_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = chembl_query_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with chembl_query_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = chembl_query_microservice.DefaultApi(api_client)
    ch_embl_molecule_request = chembl_query_microservice.ChEMBLMoleculeRequest() # ChEMBLMoleculeRequest | 

    try:
        # Query
        api_response = api_instance.query_query_chembl_post(ch_embl_molecule_request)
        print("The response of DefaultApi->query_query_chembl_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->query_query_chembl_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **ch_embl_molecule_request** | [**ChEMBLMoleculeRequest**](ChEMBLMoleculeRequest.md)|  | 

### Return type

[**ChEMBLMoleculeResponse**](ChEMBLMoleculeResponse.md)

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

