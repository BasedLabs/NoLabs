# rosettafold_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**instances_running_instances_running_get**](DefaultApi.md#instances_running_instances_running_get) | **GET** /instances-running | Instances Running
[**run_folding_run_folding_post**](DefaultApi.md#run_folding_run_folding_post) | **POST** /run-folding | Run Folding


# **instances_running_instances_running_get**
> int instances_running_instances_running_get()

Instances Running

Get number of rosettafold instances running

### Example


```python
import rosettafold_microservice
from rosettafold_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = rosettafold_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with rosettafold_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = rosettafold_microservice.DefaultApi(api_client)

    try:
        # Instances Running
        api_response = api_instance.instances_running_instances_running_get()
        print("The response of DefaultApi->instances_running_instances_running_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->instances_running_instances_running_get: %s\n" % e)
```



### Parameters

This endpoint does not need any parameter.

### Return type

**int**

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

# **run_folding_run_folding_post**
> RunRosettaFoldResponse run_folding_run_folding_post(fasta=fasta, a3m=a3m)

Run Folding

Run folding on a given amino-acid sequence :param fasta UploadFile | None: Fasta file with amino acid sequence. You must specify either fasta or a3m :param a3m UploadFile | None: MSA file. You must specify either fasta or a3m

### Example


```python
import rosettafold_microservice
from rosettafold_microservice.models.run_rosetta_fold_response import RunRosettaFoldResponse
from rosettafold_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = rosettafold_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with rosettafold_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = rosettafold_microservice.DefaultApi(api_client)
    fasta = None # bytearray |  (optional)
    a3m = None # bytearray |  (optional)

    try:
        # Run Folding
        api_response = api_instance.run_folding_run_folding_post(fasta=fasta, a3m=a3m)
        print("The response of DefaultApi->run_folding_run_folding_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->run_folding_run_folding_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **fasta** | **bytearray**|  | [optional] 
 **a3m** | **bytearray**|  | [optional] 

### Return type

[**RunRosettaFoldResponse**](RunRosettaFoldResponse.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: multipart/form-data
 - **Accept**: application/json

### HTTP response details

| Status code | Description | Response headers |
|-------------|-------------|------------------|
**200** | Successful Response |  -  |
**422** | Validation Error |  -  |

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

