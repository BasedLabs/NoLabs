# reinvent_microservice.ReinventApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**delete_api_reinvent_config_id_delete**](ReinventApi.md#delete_api_reinvent_config_id_delete) | **DELETE** /api/reinvent/{config_id} | Delete
[**get_all_configs_api_reinvent_get**](ReinventApi.md#get_all_configs_api_reinvent_get) | **GET** /api/reinvent/ | Get All Configs
[**get_config_api_reinvent_reinvent_config_id_get**](ReinventApi.md#get_config_api_reinvent_reinvent_config_id_get) | **GET** /api/reinvent/reinvent/{config_id} | Get Config
[**learning_api_reinvent_config_id_start_learning_post**](ReinventApi.md#learning_api_reinvent_config_id_start_learning_post) | **POST** /api/reinvent/{config_id}/start-learning | Learning
[**logs_api_reinvent_config_id_logs_get**](ReinventApi.md#logs_api_reinvent_config_id_logs_get) | **GET** /api/reinvent/{config_id}/logs | Logs
[**params_api_reinvent_config_id_params_get**](ReinventApi.md#params_api_reinvent_config_id_params_get) | **GET** /api/reinvent/{config_id}/params | Params
[**sampling_api_reinvent_config_id_start_sampling_post**](ReinventApi.md#sampling_api_reinvent_config_id_start_sampling_post) | **POST** /api/reinvent/{config_id}/start-sampling | Sampling
[**save_params_api_reinvent_config_id_params_post**](ReinventApi.md#save_params_api_reinvent_config_id_params_post) | **POST** /api/reinvent/{config_id}/params | Save Params
[**smiles_api_reinvent_config_id_smiles_get**](ReinventApi.md#smiles_api_reinvent_config_id_smiles_get) | **GET** /api/reinvent/{config_id}/smiles | Smiles
[**stop_api_reinvent_config_id_jobs_stop_post**](ReinventApi.md#stop_api_reinvent_config_id_jobs_stop_post) | **POST** /api/reinvent/{config_id}/jobs/stop | Stop


# **delete_api_reinvent_config_id_delete**
> object delete_api_reinvent_config_id_delete(config_id)

Delete

Delete configuration.

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
    api_instance = reinvent_microservice.ReinventApi(api_client)
    config_id = 'config_id_example' # str | 

    try:
        # Delete
        api_response = api_instance.delete_api_reinvent_config_id_delete(config_id)
        print("The response of ReinventApi->delete_api_reinvent_config_id_delete:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ReinventApi->delete_api_reinvent_config_id_delete: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **config_id** | **str**|  | 

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

# **get_all_configs_api_reinvent_get**
> List[ConfigurationResponse] get_all_configs_api_reinvent_get()

Get All Configs

Get all configurations available.

### Example


```python
import reinvent_microservice
from reinvent_microservice.models.configuration_response import ConfigurationResponse
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
    api_instance = reinvent_microservice.ReinventApi(api_client)

    try:
        # Get All Configs
        api_response = api_instance.get_all_configs_api_reinvent_get()
        print("The response of ReinventApi->get_all_configs_api_reinvent_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ReinventApi->get_all_configs_api_reinvent_get: %s\n" % e)
```



### Parameters

This endpoint does not need any parameter.

### Return type

[**List[ConfigurationResponse]**](ConfigurationResponse.md)

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

# **get_config_api_reinvent_reinvent_config_id_get**
> ResponseGetConfigApiReinventReinventConfigIdGet get_config_api_reinvent_reinvent_config_id_get(config_id)

Get Config

Get configuration.

### Example


```python
import reinvent_microservice
from reinvent_microservice.models.response_get_config_api_reinvent_reinvent_config_id_get import ResponseGetConfigApiReinventReinventConfigIdGet
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
    api_instance = reinvent_microservice.ReinventApi(api_client)
    config_id = 'config_id_example' # str | 

    try:
        # Get Config
        api_response = api_instance.get_config_api_reinvent_reinvent_config_id_get(config_id)
        print("The response of ReinventApi->get_config_api_reinvent_reinvent_config_id_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ReinventApi->get_config_api_reinvent_reinvent_config_id_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **config_id** | **str**|  | 

### Return type

[**ResponseGetConfigApiReinventReinventConfigIdGet**](ResponseGetConfigApiReinventReinventConfigIdGet.md)

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

# **learning_api_reinvent_config_id_start_learning_post**
> object learning_api_reinvent_config_id_start_learning_post(config_id)

Learning

Start model learning.

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
    api_instance = reinvent_microservice.ReinventApi(api_client)
    config_id = 'config_id_example' # str | 

    try:
        # Learning
        api_response = api_instance.learning_api_reinvent_config_id_start_learning_post(config_id)
        print("The response of ReinventApi->learning_api_reinvent_config_id_start_learning_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ReinventApi->learning_api_reinvent_config_id_start_learning_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **config_id** | **str**|  | 

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

# **logs_api_reinvent_config_id_logs_get**
> ResponseLogsApiReinventConfigIdLogsGet logs_api_reinvent_config_id_logs_get(config_id)

Logs

Get logs of all runs jobs of given configuration.

### Example


```python
import reinvent_microservice
from reinvent_microservice.models.response_logs_api_reinvent_config_id_logs_get import ResponseLogsApiReinventConfigIdLogsGet
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
    api_instance = reinvent_microservice.ReinventApi(api_client)
    config_id = 'config_id_example' # str | 

    try:
        # Logs
        api_response = api_instance.logs_api_reinvent_config_id_logs_get(config_id)
        print("The response of ReinventApi->logs_api_reinvent_config_id_logs_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ReinventApi->logs_api_reinvent_config_id_logs_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **config_id** | **str**|  | 

### Return type

[**ResponseLogsApiReinventConfigIdLogsGet**](ResponseLogsApiReinventConfigIdLogsGet.md)

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

# **params_api_reinvent_config_id_params_get**
> ResponseParamsApiReinventConfigIdParamsGet params_api_reinvent_config_id_params_get(config_id)

Params

Get learning parameters of configuration.

### Example


```python
import reinvent_microservice
from reinvent_microservice.models.response_params_api_reinvent_config_id_params_get import ResponseParamsApiReinventConfigIdParamsGet
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
    api_instance = reinvent_microservice.ReinventApi(api_client)
    config_id = 'config_id_example' # str | 

    try:
        # Params
        api_response = api_instance.params_api_reinvent_config_id_params_get(config_id)
        print("The response of ReinventApi->params_api_reinvent_config_id_params_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ReinventApi->params_api_reinvent_config_id_params_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **config_id** | **str**|  | 

### Return type

[**ResponseParamsApiReinventConfigIdParamsGet**](ResponseParamsApiReinventConfigIdParamsGet.md)

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

# **sampling_api_reinvent_config_id_start_sampling_post**
> object sampling_api_reinvent_config_id_start_sampling_post(config_id, sampling_size_request)

Sampling

Generate new ligands based on the provided config id.

### Example


```python
import reinvent_microservice
from reinvent_microservice.models.sampling_size_request import SamplingSizeRequest
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
    api_instance = reinvent_microservice.ReinventApi(api_client)
    config_id = 'config_id_example' # str | 
    sampling_size_request = reinvent_microservice.SamplingSizeRequest() # SamplingSizeRequest | 

    try:
        # Sampling
        api_response = api_instance.sampling_api_reinvent_config_id_start_sampling_post(config_id, sampling_size_request)
        print("The response of ReinventApi->sampling_api_reinvent_config_id_start_sampling_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ReinventApi->sampling_api_reinvent_config_id_start_sampling_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **config_id** | **str**|  | 
 **sampling_size_request** | [**SamplingSizeRequest**](SamplingSizeRequest.md)|  | 

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

# **save_params_api_reinvent_config_id_params_post**
> object save_params_api_reinvent_config_id_params_post(config_id, name, center_x, center_y, center_z, size_x, size_y, size_z, pdb_file, epochs=epochs, batch_size=batch_size, minscore=minscore)

Save Params

Save parameters for reinvent reinforcement learning configuration.

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
    api_instance = reinvent_microservice.ReinventApi(api_client)
    config_id = 'config_id_example' # str | 
    name = 'name_example' # str | 
    center_x = 3.4 # float | 
    center_y = 3.4 # float | 
    center_z = 3.4 # float | 
    size_x = 3.4 # float | 
    size_y = 3.4 # float | 
    size_z = 3.4 # float | 
    pdb_file = None # bytearray | 
    epochs = 50 # int |  (optional) (default to 50)
    batch_size = 128 # int |  (optional) (default to 128)
    minscore = 0.4 # float |  (optional) (default to 0.4)

    try:
        # Save Params
        api_response = api_instance.save_params_api_reinvent_config_id_params_post(config_id, name, center_x, center_y, center_z, size_x, size_y, size_z, pdb_file, epochs=epochs, batch_size=batch_size, minscore=minscore)
        print("The response of ReinventApi->save_params_api_reinvent_config_id_params_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ReinventApi->save_params_api_reinvent_config_id_params_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **config_id** | **str**|  | 
 **name** | **str**|  | 
 **center_x** | **float**|  | 
 **center_y** | **float**|  | 
 **center_z** | **float**|  | 
 **size_x** | **float**|  | 
 **size_y** | **float**|  | 
 **size_z** | **float**|  | 
 **pdb_file** | **bytearray**|  | 
 **epochs** | **int**|  | [optional] [default to 50]
 **batch_size** | **int**|  | [optional] [default to 128]
 **minscore** | **float**|  | [optional] [default to 0.4]

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

# **smiles_api_reinvent_config_id_smiles_get**
> SmilesResponse smiles_api_reinvent_config_id_smiles_get(config_id)

Smiles

Get generated smiles after sampling or RL.

### Example


```python
import reinvent_microservice
from reinvent_microservice.models.smiles_response import SmilesResponse
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
    api_instance = reinvent_microservice.ReinventApi(api_client)
    config_id = 'config_id_example' # str | 

    try:
        # Smiles
        api_response = api_instance.smiles_api_reinvent_config_id_smiles_get(config_id)
        print("The response of ReinventApi->smiles_api_reinvent_config_id_smiles_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ReinventApi->smiles_api_reinvent_config_id_smiles_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **config_id** | **str**|  | 

### Return type

[**SmilesResponse**](SmilesResponse.md)

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

# **stop_api_reinvent_config_id_jobs_stop_post**
> object stop_api_reinvent_config_id_jobs_stop_post(config_id)

Stop

Stop current job.

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
    api_instance = reinvent_microservice.ReinventApi(api_client)
    config_id = 'config_id_example' # str | 

    try:
        # Stop
        api_response = api_instance.stop_api_reinvent_config_id_jobs_stop_post(config_id)
        print("The response of ReinventApi->stop_api_reinvent_config_id_jobs_stop_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ReinventApi->stop_api_reinvent_config_id_jobs_stop_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **config_id** | **str**|  | 

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

