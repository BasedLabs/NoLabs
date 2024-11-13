import os
from transformers import EsmForProteinFolding
EsmForProteinFolding.from_pretrained('facebook/esmfold_v1', low_cpu_mem_usage=True, cache_dir=os.environ["ESMFOLD_WEIGHTS_LOCATION"])