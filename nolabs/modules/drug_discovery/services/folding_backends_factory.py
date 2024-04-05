import tempfile

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.modules.drug_discovery.services.folding_methods import FoldingMethods
from nolabs.infrastructure.settings import Settings


def run_folding(fasta_content: str, settings: Settings, method: FoldingMethods) -> str:
    FoldingMethods.ensure_contains(method)

    def run_esmfold_light():
        import esmfold_light_microservice as microservice
        from esmfold_light_microservice import ApiClient, DefaultApi, Configuration

        configuration = Configuration(
            host=settings.esmfold_light_host
        )

        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            if len(fasta_content) > 400:
                raise NoLabsException(messages=["Light folding does not support sequences longer than 400. Please use "
                                                "other folding backend"],
                                      error_code=ErrorCodes.drug_discovery_folding_error)
            request = microservice.RunEsmFoldPredictionRequest(protein_sequence=fasta_content)
            return api_instance.predict_through_api_run_folding_post(
                run_esm_fold_prediction_request=request).pdb_content

    def run_esmfold():
        import esmfold_microservice as microservice
        from esmfold_microservice import ApiClient, DefaultApi, Configuration

        configuration = Configuration(
            host=settings.esmfold_host,
        )

        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            request = microservice.RunEsmFoldPredictionRequest(protein_sequence=fasta_content)
            return api_instance.predict_run_folding_post(run_esm_fold_prediction_request=request).pdb_content

    def run_rosetta():
        from rosettafold_microservice import ApiClient, DefaultApi, Configuration

        configuration = Configuration(
            host=settings.rosettafold_host
        )
        with ApiClient(configuration=configuration) as client:
            api_instance = DefaultApi(client)
            open('tmp.fasta', 'w').write(fasta_content)
            response = api_instance.run_folding_run_folding_post(
                fasta='tmp.fasta',
                a3m=None
            )

            return response.pdb_content.anyof_schema_1_validator

    if method == FoldingMethods.esmfold_light:
        return run_esmfold_light()

    if method == FoldingMethods.esmfold:
        return run_esmfold()

    if method == FoldingMethods.rosettafold:
        return run_rosetta()




