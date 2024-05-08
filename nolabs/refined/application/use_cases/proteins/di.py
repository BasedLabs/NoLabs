from nolabs.refined.application.use_cases.proteins.use_cases import SearchProteinsFeature, GetProteinFeature, DeleteProteinFeature, UploadProteinFeature


class ProteinsControllerDependencies:
    @staticmethod
    def search_proteins() -> SearchProteinsFeature:
        return SearchProteinsFeature()

    @staticmethod
    def get_protein() -> GetProteinFeature:
        return GetProteinFeature()

    @staticmethod
    def upload_protein() -> UploadProteinFeature:
        return UploadProteinFeature()

    @staticmethod
    def delete_protein() -> DeleteProteinFeature:
        return DeleteProteinFeature()