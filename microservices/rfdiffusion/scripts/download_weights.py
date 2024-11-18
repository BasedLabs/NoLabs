import os
import urllib.request
from pathlib import Path


def download_file(url, file_name):
    if os.path.exists(file_name):
        return

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

urls = [
"http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt",
"http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt",
"http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt",
"http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt",
"http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt",
"http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt",
"http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt"
]

for url in urls:
    download_model_weights(os.environ["RFDIFFUSION_WEIGHTS_LOCATION"], url)
