# sc_gpt_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**embed_embed_post**](DefaultApi.md#embed_embed_post) | **POST** /embed | Embed
[**reference_classify_reference_classify_post**](DefaultApi.md#reference_classify_reference_classify_post) | **POST** /reference_classify | Reference Classify


# **embed_embed_post**
> EmbedResponse embed_embed_post(embed_request)

Embed

### Example


```python
import sc_gpt_microservice
from sc_gpt_microservice.models.embed_request import EmbedRequest
from sc_gpt_microservice.models.embed_response import EmbedResponse
from sc_gpt_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = sc_gpt_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with sc_gpt_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = sc_gpt_microservice.DefaultApi(api_client)
    embed_request = sc_gpt_microservice.EmbedRequest() # EmbedRequest | 

    try:
        # Embed
        api_response = api_instance.embed_embed_post(embed_request)
        print("The response of DefaultApi->embed_embed_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->embed_embed_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **embed_request** | [**EmbedRequest**](EmbedRequest.md)|  | 

### Return type

[**EmbedResponse**](EmbedResponse.md)

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

# **reference_classify_reference_classify_post**
> ReferenceMappingResponse reference_classify_reference_classify_post(reference_mapping_request)

Reference Classify

### Example


```python
import sc_gpt_microservice
from sc_gpt_microservice.models.reference_mapping_request import ReferenceMappingRequest
from sc_gpt_microservice.models.reference_mapping_response import ReferenceMappingResponse
from sc_gpt_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = sc_gpt_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with sc_gpt_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = sc_gpt_microservice.DefaultApi(api_client)
    reference_mapping_request = sc_gpt_microservice.ReferenceMappingRequest() # ReferenceMappingRequest | 

    try:
        # Reference Classify
        api_response = api_instance.reference_classify_reference_classify_post(reference_mapping_request)
        print("The response of DefaultApi->reference_classify_reference_classify_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->reference_classify_reference_classify_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **reference_mapping_request** | [**ReferenceMappingRequest**](ReferenceMappingRequest.md)|  | 

### Return type

[**ReferenceMappingResponse**](ReferenceMappingResponse.md)

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

