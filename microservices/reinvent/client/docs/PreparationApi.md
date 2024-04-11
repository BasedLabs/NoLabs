# reinvent_microservice.PreparationApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**prepare_target_api_preparation_prepare_target_post**](PreparationApi.md#prepare_target_api_preparation_prepare_target_post) | **POST** /api/preparation/prepare-target | Prepare Target


# **prepare_target_api_preparation_prepare_target_post**
> object prepare_target_api_preparation_prepare_target_post(pdb_content)

Prepare Target

Prepare .pdbqt file from pdb file.

### Example


```python
import reinvent_microservice
from reinvent_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = reinvent_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with reinvent_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = reinvent_microservice.PreparationApi(api_client)
    pdb_content = None # bytearray | 

    try:
        # Prepare Target
        api_response = api_instance.prepare_target_api_preparation_prepare_target_post(pdb_content)
        print("The response of PreparationApi->prepare_target_api_preparation_prepare_target_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling PreparationApi->prepare_target_api_preparation_prepare_target_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **pdb_content** | **bytearray**|  | 

### Return type

**object**

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

