# blast_query_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**blast_blast_post**](DefaultApi.md#blast_blast_post) | **POST** /blast | Blast
[**get_running_jobs_jobs_running_get**](DefaultApi.md#get_running_jobs_jobs_running_get) | **GET** /jobs/running | Get Running Jobs
[**is_job_running_job_job_id_is_running_get**](DefaultApi.md#is_job_running_job_job_id_is_running_get) | **GET** /job/{job_id}/is-running | Is Job Running


# **blast_blast_post**
> object blast_blast_post(sequence_query)

Blast

Perform a BLAST query with the specified type.  Args:     query (SequenceQuery): The query parameters including the sequence and BLAST type.  Returns:     Dict[str, Any]: The result of the BLAST query.

### Example


```python
import blast_query_microservice
from blast_query_microservice.models.sequence_query import SequenceQuery
from blast_query_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = blast_query_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with blast_query_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = blast_query_microservice.DefaultApi(api_client)
    sequence_query = blast_query_microservice.SequenceQuery() # SequenceQuery | 

    try:
        # Blast
        api_response = api_instance.blast_blast_post(sequence_query)
        print("The response of DefaultApi->blast_blast_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->blast_blast_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **sequence_query** | [**SequenceQuery**](SequenceQuery.md)|  | 

### Return type

**object**

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

# **get_running_jobs_jobs_running_get**
> object get_running_jobs_jobs_running_get()

Get Running Jobs

### Example


```python
import blast_query_microservice
from blast_query_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = blast_query_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with blast_query_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = blast_query_microservice.DefaultApi(api_client)

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
import blast_query_microservice
from blast_query_microservice.models.is_job_running_response import IsJobRunningResponse
from blast_query_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = blast_query_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with blast_query_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = blast_query_microservice.DefaultApi(api_client)
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

