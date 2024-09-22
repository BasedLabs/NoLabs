import base64

from rdkit import Chem
from rdkit.Chem import Draw


def generate_png_from_smiles(smiles: str | bytes) -> bytes:
    if isinstance(smiles, bytes):
        smiles = smiles.decode("utf-8")

    molecule = Chem.MolFromSmiles(smiles)
    # Draw the molecule to a PNG file
    image = Draw.MolToImage(molecule)

    from io import BytesIO

    buffer = BytesIO()
    image.save(buffer, format="PNG")
    return buffer.getvalue()


def image_file_to_base64(path):
    """
    Reads an image from the specified path and converts it to a Base64-encoded string.

    Args:
        path (str): The path to the image file.

    Returns:
        str: The Base64-encoded image data.
    """
    with open(path, "rb") as image_file:
        # Read the image data
        image_data = image_file.read()
        # Encode the data as Base64
        return base64.b64encode(image_data).decode("utf-8")
