import os
import urllib.request
from pathlib import Path

from dotenv import load_dotenv

load_dotenv(".env")

from log import logger
from settings import settings


def download_file(url, file_name):
    with urllib.request.urlopen(url) as response:
        total_size = int(response.headers['Content-Length'])

        with open(file_name, 'wb') as file:
            chunk_size = 1024 * 1024  # 1 MB
            downloaded_size = 0

            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                file.write(chunk)
                downloaded_size += len(chunk)

                logger.info(f"Downloaded: {downloaded_size / (1024 * 1024):.2f} MB / {total_size / (1024 * 1024):.2f} MB")

    logger.info(f"Download complete: {file_name}")


def download_model_weights():
    logger.info(
        "Downloading model weights", extra={"model_weight_url": settings.model1_url}
    )

    path1 = Path(settings.weights_path, settings.model1_url.split("/")[-1])
    download_file(
        settings.model1_url, str(path1)
    )

    logger.info(
        "Downloading model weights", extra={"model_weight_url": settings.model2_url}
    )

    path2 = Path(settings.weights_path, settings.model2_url.split("/")[-1])
    download_file(
        settings.model2_url, str(path2)
    )

    logger.info("Finished")

download_model_weights()