from nolabs.refined.infrastructure.settings import Settings


class InfrastructureDependencies:
    @staticmethod
    def settings() -> Settings:
        return Settings.load()
