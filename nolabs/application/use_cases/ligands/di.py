from nolabs.application.use_cases.ligands.use_cases import (SearchLigandsContentFeature,
                                                            SearchLigandsMetadataFeature,
                                                            GetLigandFeature,
                                                            DeleteLigandFeature,
                                                            UploadLigandFeature,
                                                            UpdateLigandFeature)


class LigandsControllerDependencies:
    @staticmethod
    def search_ligands_content() -> SearchLigandsContentFeature:
        return SearchLigandsContentFeature()

    @staticmethod
    def search_ligands_metadata() -> SearchLigandsMetadataFeature:
        return SearchLigandsMetadataFeature()

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
