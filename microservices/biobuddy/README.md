# BioBuddy Microservice Testing Instructions

Welcome to the BioBuddy microservice. Follow these instructions to set up and test the API.

## Prerequisites

- Python 3.11
- Make sure your openai api key is set up https://platform.openai.com/docs/quickstart/step-2-set-up-your-api-key

## Setup Instructions

### 1. Navigate to the BioBuddy Directory
Open your terminal and change to the BioBuddy directory:

```bash
cd nolabs/microservices/biobuddy
```

### 2. Create a Python Virtual Environment
Ensure you are using Python 3.11 for the virtual environment:

```bash
python -m venv venv
```

Activate the virtual environment:

For Windows:
```bash
.\venv\Scripts\activate
```

For Unix or MacOS:
```bash
source venv/bin/activate
```

### 3. Install Required Packages
Install all dependencies using pip:

```bash
pip install -r requirements.txt
```

### 4. Run the Application
Start the BioBuddy microservice using uvicorn:

```bash
uvicorn biobuddy.api:app --host 127.0.0.1 --port 5738
```

### 5. Access the API Documentation
With the server running, open your web browser and navigate to: http://127.0.0.1:5738/docs

Here you can view the Swagger documentation and interact with the API endpoints.

### 6. Restarting the Service
If you make changes to the BioBuddy service, you can restart it by stopping the running uvicorn server. Press Ctrl+C in the terminal to stop it, then restart the server with the same uvicorn command:

```bash
uvicorn biobuddy.api:app --host 127.0.0.1 --port 5738
```