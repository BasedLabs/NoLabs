import os
import urllib.request
from pathlib import Path


def download_file(url, file_name):
    with urllib.request.urlopen(url) as response:
        total_size = int(response.headers["Content-Length"])

        with open(file_name, "wb") as file:
            chunk_size = 1024 * 1024  # 1 MB
            downloaded_size = 0

            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                file.write(chunk)
                downloaded_size += len(chunk)

                print(
                    f"Downloaded: {downloaded_size / (1024 * 1024):.2f} MB / {total_size / (1024 * 1024):.2f} MB"
                )

    print(f"Download complete: {file_name}")


def download_model_weights(dir: str, model1_url: str, model2_url: str):
    print(f"Downloading model weights, model_weight_url: {model1_url}")

    path1 = Path(dir, model1_url.split("/")[-1])
    download_file(model1_url, str(path1))

    print(f"Downloading model weights, model_weight_url: {model2_url}")

    path2 = Path(dir, model2_url.split("/")[-1])
    download_file(model2_url, str(path2))

    print("Finished")


download_model_weights(os.environ[""])
