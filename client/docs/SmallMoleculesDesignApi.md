# nolabs_microservice.SmallMoleculesDesignApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**change_experiment_name_api_v1_small_molecules_design_experiment_name_post**](SmallMoleculesDesignApi.md#change_experiment_name_api_v1_small_molecules_design_experiment_name_post) | **POST** /api/v1/small-molecules-design/experiment/name | Change Experiment Name
[**create_experiment_api_v1_small_molecules_design_experiment_create_post**](SmallMoleculesDesignApi.md#create_experiment_api_v1_small_molecules_design_experiment_create_post) | **POST** /api/v1/small-molecules-design/experiment/create | Create Experiment
[**delete_api_v1_small_molecules_design_experiment_experiment_id_delete**](SmallMoleculesDesignApi.md#delete_api_v1_small_molecules_design_experiment_experiment_id_delete) | **DELETE** /api/v1/small-molecules-design/experiment/{experiment_id} | Delete
[**experiments_api_v1_small_molecules_design_experiments_get**](SmallMoleculesDesignApi.md#experiments_api_v1_small_molecules_design_experiments_get) | **GET** /api/v1/small-molecules-design/experiments | Experiments
[**get_experiment_api_v1_small_molecules_design_experiment_experiment_id_get**](SmallMoleculesDesignApi.md#get_experiment_api_v1_small_molecules_design_experiment_experiment_id_get) | **GET** /api/v1/small-molecules-design/experiment/{experiment_id} | Get Experiment
[**learning_api_v1_small_molecules_design_experiment_experiment_id_learning_post**](SmallMoleculesDesignApi.md#learning_api_v1_small_molecules_design_experiment_experiment_id_learning_post) | **POST** /api/v1/small-molecules-design/experiment/{experiment_id}/learning | Learning
[**logs_api_v1_small_molecules_design_experiment_experiment_id_logs_get**](SmallMoleculesDesignApi.md#logs_api_v1_small_molecules_design_experiment_experiment_id_logs_get) | **GET** /api/v1/small-molecules-design/experiment/{experiment_id}/logs | Logs
[**sampling_api_v1_small_molecules_design_experiment_experiment_id_sampling_post**](SmallMoleculesDesignApi.md#sampling_api_v1_small_molecules_design_experiment_experiment_id_sampling_post) | **POST** /api/v1/small-molecules-design/experiment/{experiment_id}/sampling | Sampling
[**save_properties_api_v1_small_molecules_design_experiment_experiment_id_props_post**](SmallMoleculesDesignApi.md#save_properties_api_v1_small_molecules_design_experiment_experiment_id_props_post) | **POST** /api/v1/small-molecules-design/experiment/{experiment_id}/props | Save Properties
[**smiles_api_v1_small_molecules_design_experiment_experiment_id_smiles_get**](SmallMoleculesDesignApi.md#smiles_api_v1_small_molecules_design_experiment_experiment_id_smiles_get) | **GET** /api/v1/small-molecules-design/experiment/{experiment_id}/smiles | Smiles
[**status_api_v1_small_molecules_design_experiment_experiment_id_status_get**](SmallMoleculesDesignApi.md#status_api_v1_small_molecules_design_experiment_experiment_id_status_get) | **GET** /api/v1/small-molecules-design/experiment/{experiment_id}/status | Status
[**stop_api_v1_small_molecules_design_experiment_experiment_id_stop_post**](SmallMoleculesDesignApi.md#stop_api_v1_small_molecules_design_experiment_experiment_id_stop_post) | **POST** /api/v1/small-molecules-design/experiment/{experiment_id}/stop | Stop


# **change_experiment_name_api_v1_small_molecules_design_experiment_name_post**
> object change_experiment_name_api_v1_small_molecules_design_experiment_name_post(change_experiment_name_request)

Change Experiment Name

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.change_experiment_name_request import ChangeExperimentNameRequest
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
    api_instance = nolabs_microservice.SmallMoleculesDesignApi(api_client)
    change_experiment_name_request = nolabs_microservice.ChangeExperimentNameRequest() # ChangeExperimentNameRequest | 

    try:
        # Change Experiment Name
        api_response = api_instance.change_experiment_name_api_v1_small_molecules_design_experiment_name_post(change_experiment_name_request)
        print("The response of SmallMoleculesDesignApi->change_experiment_name_api_v1_small_molecules_design_experiment_name_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling SmallMoleculesDesignApi->change_experiment_name_api_v1_small_molecules_design_experiment_name_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **change_experiment_name_request** | [**ChangeExperimentNameRequest**](ChangeExperimentNameRequest.md)|  | 

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

# **create_experiment_api_v1_small_molecules_design_experiment_create_post**
> ExperimentMetadataResponse create_experiment_api_v1_small_molecules_design_experiment_create_post()

Create Experiment

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.experiment_metadata_response import ExperimentMetadataResponse
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
    api_instance = nolabs_microservice.SmallMoleculesDesignApi(api_client)

    try:
        # Create Experiment
        api_response = api_instance.create_experiment_api_v1_small_molecules_design_experiment_create_post()
        print("The response of SmallMoleculesDesignApi->create_experiment_api_v1_small_molecules_design_experiment_create_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling SmallMoleculesDesignApi->create_experiment_api_v1_small_molecules_design_experiment_create_post: %s\n" % e)
```



### Parameters

This endpoint does not need any parameter.

### Return type

[**ExperimentMetadataResponse**](ExperimentMetadataResponse.md)

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

# **delete_api_v1_small_molecules_design_experiment_experiment_id_delete**
> object delete_api_v1_small_molecules_design_experiment_experiment_id_delete(experiment_id)

Delete

### Example


```python
import time
import os
import nolabs_microservice
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
    api_instance = nolabs_microservice.SmallMoleculesDesignApi(api_client)
    experiment_id = None # object | 

    try:
        # Delete
        api_response = api_instance.delete_api_v1_small_molecules_design_experiment_experiment_id_delete(experiment_id)
        print("The response of SmallMoleculesDesignApi->delete_api_v1_small_molecules_design_experiment_experiment_id_delete:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling SmallMoleculesDesignApi->delete_api_v1_small_molecules_design_experiment_experiment_id_delete: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | [**object**](.md)|  | 

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

# **experiments_api_v1_small_molecules_design_experiments_get**
> object experiments_api_v1_small_molecules_design_experiments_get()

Experiments

### Example


```python
import time
import os
import nolabs_microservice
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
    api_instance = nolabs_microservice.SmallMoleculesDesignApi(api_client)

    try:
        # Experiments
        api_response = api_instance.experiments_api_v1_small_molecules_design_experiments_get()
        print("The response of SmallMoleculesDesignApi->experiments_api_v1_small_molecules_design_experiments_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling SmallMoleculesDesignApi->experiments_api_v1_small_molecules_design_experiments_get: %s\n" % e)
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

# **get_experiment_api_v1_small_molecules_design_experiment_experiment_id_get**
> NolabsApiModelsSmallMoleculesDesignGetExperimentResponse get_experiment_api_v1_small_molecules_design_experiment_experiment_id_get(experiment_id)

Get Experiment

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.nolabs_api_models_small_molecules_design_get_experiment_response import NolabsApiModelsSmallMoleculesDesignGetExperimentResponse
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
    api_instance = nolabs_microservice.SmallMoleculesDesignApi(api_client)
    experiment_id = None # object | 

    try:
        # Get Experiment
        api_response = api_instance.get_experiment_api_v1_small_molecules_design_experiment_experiment_id_get(experiment_id)
        print("The response of SmallMoleculesDesignApi->get_experiment_api_v1_small_molecules_design_experiment_experiment_id_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling SmallMoleculesDesignApi->get_experiment_api_v1_small_molecules_design_experiment_experiment_id_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | [**object**](.md)|  | 

### Return type

[**NolabsApiModelsSmallMoleculesDesignGetExperimentResponse**](NolabsApiModelsSmallMoleculesDesignGetExperimentResponse.md)

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

# **learning_api_v1_small_molecules_design_experiment_experiment_id_learning_post**
> object learning_api_v1_small_molecules_design_experiment_experiment_id_learning_post(experiment_id)

Learning

### Example


```python
import time
import os
import nolabs_microservice
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
    api_instance = nolabs_microservice.SmallMoleculesDesignApi(api_client)
    experiment_id = None # object | 

    try:
        # Learning
        api_response = api_instance.learning_api_v1_small_molecules_design_experiment_experiment_id_learning_post(experiment_id)
        print("The response of SmallMoleculesDesignApi->learning_api_v1_small_molecules_design_experiment_experiment_id_learning_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling SmallMoleculesDesignApi->learning_api_v1_small_molecules_design_experiment_experiment_id_learning_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | [**object**](.md)|  | 

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

# **logs_api_v1_small_molecules_design_experiment_experiment_id_logs_get**
> LogsResponse logs_api_v1_small_molecules_design_experiment_experiment_id_logs_get(experiment_id)

Logs

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.logs_response import LogsResponse
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
    api_instance = nolabs_microservice.SmallMoleculesDesignApi(api_client)
    experiment_id = None # object | 

    try:
        # Logs
        api_response = api_instance.logs_api_v1_small_molecules_design_experiment_experiment_id_logs_get(experiment_id)
        print("The response of SmallMoleculesDesignApi->logs_api_v1_small_molecules_design_experiment_experiment_id_logs_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling SmallMoleculesDesignApi->logs_api_v1_small_molecules_design_experiment_experiment_id_logs_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | [**object**](.md)|  | 

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

# **sampling_api_v1_small_molecules_design_experiment_experiment_id_sampling_post**
> object sampling_api_v1_small_molecules_design_experiment_experiment_id_sampling_post(experiment_id, sampling_size_request)

Sampling

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.sampling_size_request import SamplingSizeRequest
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
    api_instance = nolabs_microservice.SmallMoleculesDesignApi(api_client)
    experiment_id = None # object | 
    sampling_size_request = nolabs_microservice.SamplingSizeRequest() # SamplingSizeRequest | 

    try:
        # Sampling
        api_response = api_instance.sampling_api_v1_small_molecules_design_experiment_experiment_id_sampling_post(experiment_id, sampling_size_request)
        print("The response of SmallMoleculesDesignApi->sampling_api_v1_small_molecules_design_experiment_experiment_id_sampling_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling SmallMoleculesDesignApi->sampling_api_v1_small_molecules_design_experiment_experiment_id_sampling_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | [**object**](.md)|  | 
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

# **save_properties_api_v1_small_molecules_design_experiment_experiment_id_props_post**
> object save_properties_api_v1_small_molecules_design_experiment_experiment_id_props_post(experiment_id, center_x, center_y, center_z, size_x, size_y, size_z, pdb_file, epochs=epochs, batch_size=batch_size, minscore=minscore)

Save Properties

### Example


```python
import time
import os
import nolabs_microservice
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
    api_instance = nolabs_microservice.SmallMoleculesDesignApi(api_client)
    experiment_id = None # object | 
    center_x = None # object | 
    center_y = None # object | 
    center_z = None # object | 
    size_x = None # object | 
    size_y = None # object | 
    size_z = None # object | 
    pdb_file = None # object | 
    epochs = None # object |  (optional)
    batch_size = None # object |  (optional)
    minscore = None # object |  (optional)

    try:
        # Save Properties
        api_response = api_instance.save_properties_api_v1_small_molecules_design_experiment_experiment_id_props_post(experiment_id, center_x, center_y, center_z, size_x, size_y, size_z, pdb_file, epochs=epochs, batch_size=batch_size, minscore=minscore)
        print("The response of SmallMoleculesDesignApi->save_properties_api_v1_small_molecules_design_experiment_experiment_id_props_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling SmallMoleculesDesignApi->save_properties_api_v1_small_molecules_design_experiment_experiment_id_props_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | [**object**](.md)|  | 
 **center_x** | [**object**](.md)|  | 
 **center_y** | [**object**](.md)|  | 
 **center_z** | [**object**](.md)|  | 
 **size_x** | [**object**](.md)|  | 
 **size_y** | [**object**](.md)|  | 
 **size_z** | [**object**](.md)|  | 
 **pdb_file** | [**object**](object.md)|  | 
 **epochs** | [**object**](.md)|  | [optional] 
 **batch_size** | [**object**](.md)|  | [optional] 
 **minscore** | [**object**](.md)|  | [optional] 

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

# **smiles_api_v1_small_molecules_design_experiment_experiment_id_smiles_get**
> object smiles_api_v1_small_molecules_design_experiment_experiment_id_smiles_get(experiment_id)

Smiles

### Example


```python
import time
import os
import nolabs_microservice
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
    api_instance = nolabs_microservice.SmallMoleculesDesignApi(api_client)
    experiment_id = None # object | 

    try:
        # Smiles
        api_response = api_instance.smiles_api_v1_small_molecules_design_experiment_experiment_id_smiles_get(experiment_id)
        print("The response of SmallMoleculesDesignApi->smiles_api_v1_small_molecules_design_experiment_experiment_id_smiles_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling SmallMoleculesDesignApi->smiles_api_v1_small_molecules_design_experiment_experiment_id_smiles_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | [**object**](.md)|  | 

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

# **status_api_v1_small_molecules_design_experiment_experiment_id_status_get**
> GetExperimentStatusResponse status_api_v1_small_molecules_design_experiment_experiment_id_status_get(experiment_id)

Status

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.get_experiment_status_response import GetExperimentStatusResponse
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
    api_instance = nolabs_microservice.SmallMoleculesDesignApi(api_client)
    experiment_id = None # object | 

    try:
        # Status
        api_response = api_instance.status_api_v1_small_molecules_design_experiment_experiment_id_status_get(experiment_id)
        print("The response of SmallMoleculesDesignApi->status_api_v1_small_molecules_design_experiment_experiment_id_status_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling SmallMoleculesDesignApi->status_api_v1_small_molecules_design_experiment_experiment_id_status_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | [**object**](.md)|  | 

### Return type

[**GetExperimentStatusResponse**](GetExperimentStatusResponse.md)

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

# **stop_api_v1_small_molecules_design_experiment_experiment_id_stop_post**
> object stop_api_v1_small_molecules_design_experiment_experiment_id_stop_post(experiment_id)

Stop

### Example


```python
import time
import os
import nolabs_microservice
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
    api_instance = nolabs_microservice.SmallMoleculesDesignApi(api_client)
    experiment_id = None # object | 

    try:
        # Stop
        api_response = api_instance.stop_api_v1_small_molecules_design_experiment_experiment_id_stop_post(experiment_id)
        print("The response of SmallMoleculesDesignApi->stop_api_v1_small_molecules_design_experiment_experiment_id_stop_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling SmallMoleculesDesignApi->stop_api_v1_small_molecules_design_experiment_experiment_id_stop_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | [**object**](.md)|  | 

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

