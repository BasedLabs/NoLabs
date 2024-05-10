from nolabs.refined.application.use_cases.ligands.use_cases import SearchLigandsFeature, GetLigandFeature, \
    DeleteLigandFeature, UploadLigandFeature, UpdateLigandFeature


class LigandsControllerDependencies:
    @staticmethod
    def search_ligands() -> SearchLigandsFeature:
        return SearchLigandsFeature()

    @staticmethod
    def get_ligand() -> GetLigandFeature:
        return GetLigandFeature()

    @staticmethod
    def upload_ligand() -> UploadLigandFeature:
        return UploadLigandFeature()

    @staticmethod
    def delete_ligand() -> DeleteLigandFeature:
        return DeleteLigandFeature()

    @staticmethod
    def update_ligand() -> UpdateLigandFeature:
        return UpdateLigandFeature()