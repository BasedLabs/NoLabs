import os

import gdown


def download_model(checkpoint_path):
    if not os.path.exists(checkpoint_path):
        os.makedirs(checkpoint_path, exist_ok=True)
        gdown.download_folder(
            "https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y",
            output=checkpoint_path,
        )


if __name__ == "__main__":
    download_model("checkpoints/scGPT_human")
