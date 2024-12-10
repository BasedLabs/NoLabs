@echo off
REM Load variables from .env if it exists
if exist .env (
    for /f "usebackq tokens=1,* delims==" %%i in (".env") do (
        set %%i=%%j
    )
)

set TARGET=%1

if "%TARGET%"=="" (
    echo Usage: build.bat [target]
    echo Available targets:
    echo   flower
    echo   install-openapi-generator
    echo   generate-client
    echo   install-mock-server
    echo   start-mock-server
    echo   gen-envs
    echo   download-diffdock-weights
    echo   download-esmfold-weights
    echo   download-rfdiffusion-weights
    echo   download-arxiv-abstracts-db
    exit /b 1
)

if "%TARGET%"=="flower" (
    celery --broker=%REDIS_URL% flower --port=5555
    goto :EOF
)

if "%TARGET%"=="install-openapi-generator" (
    npm install -g openapi-typescript-codegen
    goto :EOF
)

if "%TARGET%"=="generate-client" (
    openapi --input "http://127.0.0.1:8000/openapi.json" --output frontend/src/api/client --client axios
    goto :EOF
)

if "%TARGET%"=="install-mock-server" (
    npm install -g @stoplight/prism-cli
    goto :EOF
)

if "%TARGET%"=="start-mock-server" (
    prism mock http://127.0.0.1:%UVICORN_HOST%/openapi.json
    goto :EOF
)

if "%TARGET%"=="gen-envs" (
    python3 scripts/gen_envs.py
    goto :EOF
)

if "%TARGET%"=="download-diffdock-weights" (
    echo Downloading diffdock model weights...
    if not exist "%DIFFDOCK_WEIGHTS_LOCATION%" mkdir "%DIFFDOCK_WEIGHTS_LOCATION%"
    python3 microservices/diffdock/scripts/download_weights.py
    echo Download complete!
    goto :EOF
)

if "%TARGET%"=="download-esmfold-weights" (
    echo Downloading esmfold model weights...
    if not exist "%ESMFOLD_WEIGHTS_LOCATION%" mkdir "%ESMFOLD_WEIGHTS_LOCATION%"
    pip3 install transformers[torch] --verbose
    python3 microservices/esmfold/scripts/download_weights.py
    echo Download complete!
    goto :EOF
)

if "%TARGET%"=="download-rfdiffusion-weights" (
    echo Downloading rfdiffusion model weights...
    if not exist "%RFDIFFUSION_WEIGHTS_LOCATION%" mkdir "%RFDIFFUSION_WEIGHTS_LOCATION%"
    python3 microservices/rfdiffusion/scripts/download_weights.py
    echo Download complete!
    goto :EOF
)

if "%TARGET%"=="download-arxiv-abstracts-db" (
    echo Downloading arxiv abstracts
    where unzip >nul 2>nul
    if errorlevel 1 (
        echo You must install unzip. For example:
        echo   choco install unzip (if using Chocolatey)
        goto :EOF
    )
    if not exist "%ARXIV_ABSTRACTS_DB%" mkdir "%ARXIV_ABSTRACTS_DB%"
    curl -L -o "%ARXIV_ABSTRACTS_DB%\chroma_db.zip" https://www.kaggle.com/api/v1/datasets/download/timurishmuratov/nolabs-arxiv-abstract-vector-db
    powershell -Command "Expand-Archive '%ARXIV_ABSTRACTS_DB%\chroma_db.zip' '%ARXIV_ABSTRACTS_DB%'"
    del "%ARXIV_ABSTRACTS_DB%\chroma_db.zip"
    echo Download complete!
    goto :EOF
)

echo Unknown target "%TARGET%"
exit /b 1
