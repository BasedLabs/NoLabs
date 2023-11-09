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
    def get_experiment(self, request):
        pass

    @abstractmethod
    def change_experiment_name(self, request: Request):
        pass

    @abstractmethod
    def delete_experiment(self, request: Request):
        pass
