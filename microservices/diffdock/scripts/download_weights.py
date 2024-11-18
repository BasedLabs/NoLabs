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


def download_model_weights(dir: str, model_url: str):
    print(f"Downloading model weights, model_weight_url: {model_url}")

    path1 = Path(dir, model_url.split("/")[-1])
    download_file(model_url, str(path1))

    print("Finished")


download_model_weights(os.environ["DIFFDOCK_WEIGHTS_LOCATION"], "https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt")
download_model_weights(os.environ["DIFFDOCK_WEIGHTS_LOCATION"], "https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50D-contact-regression.pt")
