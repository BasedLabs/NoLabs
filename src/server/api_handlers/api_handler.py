import uuid
from abc import abstractmethod

from flask import Request

from src.server.services.mixins import UUIDGenerator


class ApiHandler(UUIDGenerator):
    @abstractmethod
    def inference(self, request: Request) -> dict:
        pass

    @abstractmethod
    def get_experiments(self):
        pass

    @abstractmethod
    def get_predictions(self):
        """
        Pulls predictions
        """
        pass

    @abstractmethod
    def get_experiment(self, request):
        """
        Pulls all experiment metadata
        """
        pass

    @abstractmethod
    def get_experiment_progress(self, request):
        """
        Pulls all experiment metadata
        """
        pass

    @abstractmethod
    def get_instance_progress(self, request):
        """
        Pulls all experiment metadata
        """
        pass

    @abstractmethod
    def change_experiment_name(self, request: Request):
        pass

    @abstractmethod
    def delete_experiment(self, request: Request):
        pass
