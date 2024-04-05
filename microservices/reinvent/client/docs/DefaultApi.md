# reinvent_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**delete_jobs_job_id_delete**](DefaultApi.md#delete_jobs_job_id_delete) | **DELETE** /jobs/{job_id} | Delete
[**get_all_jobs_jobs_get**](DefaultApi.md#get_all_jobs_jobs_get) | **GET** /jobs | Get All Jobs
[**get_job_jobs_job_id_get**](DefaultApi.md#get_job_jobs_job_id_get) | **GET** /jobs/{job_id} | Get Job
[**logs_jobs_job_id_logs_get**](DefaultApi.md#logs_jobs_job_id_logs_get) | **GET** /jobs/{job_id}/logs | Logs
[**params_jobs_job_id_params_get**](DefaultApi.md#params_jobs_job_id_params_get) | **GET** /jobs/{job_id}/params | Params
[**prepare_binder_prepare_binder_post**](DefaultApi.md#prepare_binder_prepare_binder_post) | **POST** /prepare-binder | Prepare Binder
[**run_jobs_job_id_run_post**](DefaultApi.md#run_jobs_job_id_run_post) | **POST** /jobs/{job_id}/run | Run
[**save_params_jobs_job_id_params_post**](DefaultApi.md#save_params_jobs_job_id_params_post) | **POST** /jobs/{job_id}/params | Save Params
[**smiles_jobs_job_id_smiles_get**](DefaultApi.md#smiles_jobs_job_id_smiles_get) | **GET** /jobs/{job_id}/smiles | Smiles
[**stop_jobs_job_id_stop_post**](DefaultApi.md#stop_jobs_job_id_stop_post) | **POST** /jobs/{job_id}/stop | Stop


# **delete_jobs_job_id_delete**
> object delete_jobs_job_id_delete(job_id)

Delete

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
    api_instance = reinvent_microservice.DefaultApi(api_client)
    job_id = 'job_id_example' # str | 

    try:
        # Delete
        api_response = api_instance.delete_jobs_job_id_delete(job_id)
        print("The response of DefaultApi->delete_jobs_job_id_delete:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->delete_jobs_job_id_delete: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **job_id** | **str**|  | 

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

# **get_all_jobs_jobs_get**
> List[JobResponse] get_all_jobs_jobs_get()

Get All Jobs

### Example


```python
import reinvent_microservice
from reinvent_microservice.models.job_response import JobResponse
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
    api_instance = reinvent_microservice.DefaultApi(api_client)

    try:
        # Get All Jobs
        api_response = api_instance.get_all_jobs_jobs_get()
        print("The response of DefaultApi->get_all_jobs_jobs_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->get_all_jobs_jobs_get: %s\n" % e)
```



### Parameters

This endpoint does not need any parameter.

### Return type

[**List[JobResponse]**](JobResponse.md)

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

# **get_job_jobs_job_id_get**
> JobResponse get_job_jobs_job_id_get(job_id)

Get Job

### Example


```python
import reinvent_microservice
from reinvent_microservice.models.job_response import JobResponse
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
    api_instance = reinvent_microservice.DefaultApi(api_client)
    job_id = 'job_id_example' # str | 

    try:
        # Get Job
        api_response = api_instance.get_job_jobs_job_id_get(job_id)
        print("The response of DefaultApi->get_job_jobs_job_id_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->get_job_jobs_job_id_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **job_id** | **str**|  | 

### Return type

[**JobResponse**](JobResponse.md)

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

# **logs_jobs_job_id_logs_get**
> LogsResponse logs_jobs_job_id_logs_get(job_id)

Logs

### Example


```python
import reinvent_microservice
from reinvent_microservice.models.logs_response import LogsResponse
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
    api_instance = reinvent_microservice.DefaultApi(api_client)
    job_id = 'job_id_example' # str | 

    try:
        # Logs
        api_response = api_instance.logs_jobs_job_id_logs_get(job_id)
        print("The response of DefaultApi->logs_jobs_job_id_logs_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->logs_jobs_job_id_logs_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **job_id** | **str**|  | 

### Return type

[**LogsResponse**](LogsResponse.md)

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

# **params_jobs_job_id_params_get**
> ParamsResponse params_jobs_job_id_params_get(job_id)

Params

### Example


```python
import reinvent_microservice
from reinvent_microservice.models.params_response import ParamsResponse
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
    api_instance = reinvent_microservice.DefaultApi(api_client)
    job_id = 'job_id_example' # str | 

    try:
        # Params
        api_response = api_instance.params_jobs_job_id_params_get(job_id)
        print("The response of DefaultApi->params_jobs_job_id_params_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->params_jobs_job_id_params_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **job_id** | **str**|  | 

### Return type

[**ParamsResponse**](ParamsResponse.md)

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

# **prepare_binder_prepare_binder_post**
> object prepare_binder_prepare_binder_post(pdb_content)

Prepare Binder

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
    api_instance = reinvent_microservice.DefaultApi(api_client)
    pdb_content = None # bytearray | 

    try:
        # Prepare Binder
        api_response = api_instance.prepare_binder_prepare_binder_post(pdb_content)
        print("The response of DefaultApi->prepare_binder_prepare_binder_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->prepare_binder_prepare_binder_post: %s\n" % e)
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

# **run_jobs_job_id_run_post**
> object run_jobs_job_id_run_post(job_id)

Run

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
    api_instance = reinvent_microservice.DefaultApi(api_client)
    job_id = 'job_id_example' # str | 

    try:
        # Run
        api_response = api_instance.run_jobs_job_id_run_post(job_id)
        print("The response of DefaultApi->run_jobs_job_id_run_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->run_jobs_job_id_run_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **job_id** | **str**|  | 

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

# **save_params_jobs_job_id_params_post**
> object save_params_jobs_job_id_params_post(job_id, name, center_x, center_y, center_z, size_x, size_y, size_z, pdb_file, epochs=epochs, batch_size=batch_size, minscore=minscore)

Save Params

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
    api_instance = reinvent_microservice.DefaultApi(api_client)
    job_id = 'job_id_example' # str | 
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
        api_response = api_instance.save_params_jobs_job_id_params_post(job_id, name, center_x, center_y, center_z, size_x, size_y, size_z, pdb_file, epochs=epochs, batch_size=batch_size, minscore=minscore)
        print("The response of DefaultApi->save_params_jobs_job_id_params_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->save_params_jobs_job_id_params_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **job_id** | **str**|  | 
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

# **smiles_jobs_job_id_smiles_get**
> SmilesResponse smiles_jobs_job_id_smiles_get(job_id)

Smiles

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
    api_instance = reinvent_microservice.DefaultApi(api_client)
    job_id = 'job_id_example' # str | 

    try:
        # Smiles
        api_response = api_instance.smiles_jobs_job_id_smiles_get(job_id)
        print("The response of DefaultApi->smiles_jobs_job_id_smiles_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->smiles_jobs_job_id_smiles_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **job_id** | **str**|  | 

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

# **stop_jobs_job_id_stop_post**
> object stop_jobs_job_id_stop_post(job_id)

Stop

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
    api_instance = reinvent_microservice.DefaultApi(api_client)
    job_id = 'job_id_example' # str | 

    try:
        # Stop
        api_response = api_instance.stop_jobs_job_id_stop_post(job_id)
        print("The response of DefaultApi->stop_jobs_job_id_stop_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->stop_jobs_job_id_stop_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **job_id** | **str**|  | 

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

