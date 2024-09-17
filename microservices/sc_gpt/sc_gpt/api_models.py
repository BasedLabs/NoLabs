import dataclasses
from typing import List, Optional


@dataclasses.dataclass
class ReferenceMappingRequest:
    """
    Represents a request for reference mapping. This will embed all cells by their gene expression profiles and then use a k-nearest neighbors classifier to predict the cell type of the query dataset using similar cells in the reference dataset.

    Attributes:
        query_path (str): Path to the query dataset in h5ad format.
        reference_path (str): Path to the reference dataset in h5ad format.
        cell_type_key (Optional[str], optional): Key to extract ground truth cell type from the query dataset. Defaults to "celltype". Reference dataset should have this key. Query dataset should have this key if you want to calculate metrics. Defaults to "celltype".
        gene_col (str, optional): Column name in the query and reference datasets that contains the gene names. Defaults to "gene_name".
        batch_size (int, optional): Batch size for processing. Defaults to 64.
        device (str, optional): Device to use for computation (e.g., "cpu", "cuda"). Defaults to "cpu".
        k_neighbors (int, optional): Number of neighbors to consider for classification task. Defaults to 10.
    """

    query_path: str
    reference_path: str
    cell_type_key: Optional[str] = "celltype"
    gene_col: str = "gene_name"
    batch_size: int = 64
    device: str = "cpu"
    k_neighbors: int = 10
    calculate_metrics: bool = False


@dataclasses.dataclass
class ReferenceMappingResponse:
    predictions: List[str]
    accuracy: Optional[float]
    precision: Optional[float]
    recall: Optional[float]
    macro_f1: Optional[float]


@dataclasses.dataclass
class EmbedRequest:
    dataset_path: str
    gene_col: str = "gene_name"
    batch_size: int = 64
    device: str = "cpu"


@dataclasses.dataclass
class EmbedResponse:
    embeddings: List[List[float]]
