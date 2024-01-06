from umol.model import DrugTargetInteraction
from umol.api_models import RunUmolPredictionResponce, RunUmolPredictionRequest
from umol.loggers import Log

__all__ = ['run_umol']

def run_umol(parameters: RunUmolPredictionRequest) -> RunUmolPredictionResponce:
    try:
        model = DrugTargetInteraction()
        sdf_content, plddt_array = model.predict(
            protein_fasta_file=parameters.protein_file,
            ligand_smiles=parameters.ligand_smiles,
            msa_file=parameters.msa_file,
            binding_pocket=parameters.binding_pocket
        )
        return RunUmolPredictionResponce(sdf_contents=sdf_content, plddt_array=plddt_array)
    except Exception as e:
        Log.exception()
        return RunUmolPredictionResponce(sdf_contents=None, plddt_array=None,
                                         errors=['Unable to run umol. Internal error', str(e)])
