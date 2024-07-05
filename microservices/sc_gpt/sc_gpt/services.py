import os
import gdown
from pathlib import Path
import numpy as np
from scipy.stats import mode
import scanpy as sc
import warnings
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import confusion_matrix
import seaborn as sns
import pandas as pd
import sys
import faiss
import scgpt as scg
from anndata import AnnData

from sc_gpt.api_models import ReferenceMappingRequest, ReferenceMappingResponse, EmbedRequest, EmbedResponse
from sc_gpt.cell_emb import embed_data

def download_model(checkpoint_path="checkpoints/scGPT_human"):
    if not os.path.exists(checkpoint_path):
        os.makedirs(checkpoint_path, exist_ok=True)
        gdown.download_folder(
            "https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y",
            output=checkpoint_path,
        )

def l2_sim(a, b):
    sims = -np.linalg.norm(a - b, axis=1)
    return sims

def get_similar_vectors(vector, ref, top_k=10):
        # sims = cos_sim(vector, ref)
        sims = l2_sim(vector, ref)

        top_k_idx = np.argsort(sims)[::-1][:top_k]
        return top_k_idx, sims[top_k_idx]

def classify_cells(request: ReferenceMappingRequest, checkpoint_path="checkpoints/scGPT_human") -> ReferenceMappingResponse:
    assert os.path.exists(request.query_path), f"Query path {request.query_path} does not exist"
    assert os.path.exists(request.reference_path), f"Reference path {request.reference_path} does not exist"
    assert os.path.exists(checkpoint_path), f"Checkpoint path {checkpoint_path} does not exist"

    placeholder_cell_type = "To be predicted"

    adata = sc.read_h5ad(request.reference_path)
    test_adata = sc.read_h5ad(request.query_path)

    # due to issue similar to https://github.com/scverse/anndata/issues/1494
    # Basically I can not change adata.X without creating a new AnnData object
    adata = AnnData(X=adata.X.toarray(), obs=adata.obs, var=adata.var, obsm=adata.obsm, varm=adata.varm, obsp=adata.obsp, varp=adata.varp, layers=adata.layers)

    # if cell type key is not present, add a placeholder for query dataset
    if test_adata.obs.get(request.cell_type_key) is None:
        test_adata.obs[request.cell_type_key] = ["To be predicted"] * len(test_adata)

    model_dir = Path(checkpoint_path)

    ref_embed_adata = embed_data(
        adata,
        model_dir,
        gene_col=request.gene_col,
        batch_size=request.batch_size,
        device=request.device,
    )

    test_embed_adata = embed_data(
        test_adata,
        model_dir,
        gene_col=request.gene_col,
        batch_size=request.batch_size,
        device=request.device,
    )

    # concatenate the two datasets
    adata_concat = test_embed_adata.concatenate(ref_embed_adata, batch_key="dataset")
    # mark the reference vs. query dataset
    adata_concat.obs["is_ref"] = ["Query"] * len(test_embed_adata) + ["Reference"] * len(
        ref_embed_adata
    )
    adata_concat.obs["is_ref"] = adata_concat.obs["is_ref"].astype("category")

    adata_concat.obs[request.cell_type_key] = adata_concat.obs[request.cell_type_key].astype("category")
    adata_concat.obs[request.cell_type_key] = adata_concat.obs[request.cell_type_key].cat.add_categories(["To be predicted"])
    adata_concat.obs[request.cell_type_key][: len(test_embed_adata)] = placeholder_cell_type

    ref_cell_embeddings = ref_embed_adata.obsm["X_scGPT"]
    test_emebd = test_embed_adata.obsm["X_scGPT"]

    index = faiss.IndexFlatL2(ref_cell_embeddings.shape[1])
    index.add(ref_cell_embeddings)

    # Query dataset, k - number of closest elements (returns 2 numpy arrays)
    distances, labels = index.search(test_emebd, request.k_neighbors)

    idx_list=[i for i in range(test_emebd.shape[0])]
    preds = []
    for k in idx_list:
        idx = labels[k]
        pred = ref_embed_adata.obs[request.cell_type_key][idx].value_counts()
        preds.append(pred.index[0])

    if request.calculate_metrics:
        gt = test_adata.obs[request.cell_type_key].to_numpy()
        accuracy = accuracy_score(gt, preds)
        precision = precision_score(gt, preds, average="macro")
        recall = recall_score(gt, preds, average="macro")
        macro_f1 = f1_score(gt, preds, average="macro")
    else:
        accuracy = None
        precision = None
        recall = None
        macro_f1 = None

    response = ReferenceMappingResponse(
        predictions=preds,
        accuracy=accuracy,
        precision=precision,
        recall=recall,
        macro_f1=macro_f1,
    )

    return response

def get_embeddings(request: EmbedRequest, checkpoint_path="checkpoints/scGPT_human") -> EmbedResponse:
    assert os.path.exists(request.dataset_path), f"Dataset path {request.dataset_path} does not exist"
    assert os.path.exists(checkpoint_path), f"Checkpoint path {checkpoint_path} does not exist"

    adata = sc.read_h5ad(request.dataset_path)

    # due to issue similar to https://github.com/scverse/anndata/issues/1494
    # Basically I can not change adata.X without creating a new AnnData object
    adata = AnnData(X=adata.X.toarray(), obs=adata.obs, var=adata.var, obsm=adata.obsm, varm=adata.varm, obsp=adata.obsp, varp=adata.varp, layers=adata.layers)

    model_dir = Path(checkpoint_path)

    ref_embed_adata = embed_data(
        adata,
        model_dir,
        gene_col=request.gene_col,
        batch_size=request.batch_size,
        device=request.device,
    )

    response = EmbedResponse(embeddings=ref_embed_adata.obsm["X_scGPT"].tolist())

    return response

