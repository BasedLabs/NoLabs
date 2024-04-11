# nolabs_microservice.BiobuddyApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**check_biobuddy_enabled_api_v1_biobuddy_check_biobuddy_enabled_get**](BiobuddyApi.md#check_biobuddy_enabled_api_v1_biobuddy_check_biobuddy_enabled_get) | **GET** /api/v1/biobuddy/check_biobuddy_enabled | Check Biobuddy Enabled
[**load_conversation_api_v1_biobuddy_load_conversation_get**](BiobuddyApi.md#load_conversation_api_v1_biobuddy_load_conversation_get) | **GET** /api/v1/biobuddy/load-conversation | Load Conversation
[**send_message_api_v1_biobuddy_send_message_post**](BiobuddyApi.md#send_message_api_v1_biobuddy_send_message_post) | **POST** /api/v1/biobuddy/send-message | Send Message


# **check_biobuddy_enabled_api_v1_biobuddy_check_biobuddy_enabled_get**
> CheckBioBuddyEnabledResponse check_biobuddy_enabled_api_v1_biobuddy_check_biobuddy_enabled_get()

Check Biobuddy Enabled

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.check_bio_buddy_enabled_response import CheckBioBuddyEnabledResponse
from nolabs_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = nolabs_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with nolabs_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = nolabs_microservice.BiobuddyApi(api_client)

    try:
        # Check Biobuddy Enabled
        api_response = api_instance.check_biobuddy_enabled_api_v1_biobuddy_check_biobuddy_enabled_get()
        print("The response of BiobuddyApi->check_biobuddy_enabled_api_v1_biobuddy_check_biobuddy_enabled_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling BiobuddyApi->check_biobuddy_enabled_api_v1_biobuddy_check_biobuddy_enabled_get: %s\n" % e)
```



### Parameters

This endpoint does not need any parameter.

### Return type

[**CheckBioBuddyEnabledResponse**](CheckBioBuddyEnabledResponse.md)

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

# **load_conversation_api_v1_biobuddy_load_conversation_get**
> LoadConversationResponse load_conversation_api_v1_biobuddy_load_conversation_get(experiment_id)

Load Conversation

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.load_conversation_response import LoadConversationResponse
from nolabs_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = nolabs_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with nolabs_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = nolabs_microservice.BiobuddyApi(api_client)
    experiment_id = 'experiment_id_example' # str | 

    try:
        # Load Conversation
        api_response = api_instance.load_conversation_api_v1_biobuddy_load_conversation_get(experiment_id)
        print("The response of BiobuddyApi->load_conversation_api_v1_biobuddy_load_conversation_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling BiobuddyApi->load_conversation_api_v1_biobuddy_load_conversation_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | **str**|  | 

### Return type

[**LoadConversationResponse**](LoadConversationResponse.md)

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

# **send_message_api_v1_biobuddy_send_message_post**
> SendMessageResponse send_message_api_v1_biobuddy_send_message_post(experiment_id, message_content)

Send Message

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.send_message_response import SendMessageResponse
from nolabs_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = nolabs_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with nolabs_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = nolabs_microservice.BiobuddyApi(api_client)
    experiment_id = 'experiment_id_example' # str | 
    message_content = 'message_content_example' # str | 

    try:
        # Send Message
        api_response = api_instance.send_message_api_v1_biobuddy_send_message_post(experiment_id, message_content)
        print("The response of BiobuddyApi->send_message_api_v1_biobuddy_send_message_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling BiobuddyApi->send_message_api_v1_biobuddy_send_message_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | **str**|  | 
 **message_content** | **str**|  | 

### Return type

[**SendMessageResponse**](SendMessageResponse.md)

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

