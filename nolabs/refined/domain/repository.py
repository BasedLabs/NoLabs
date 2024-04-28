from mongoengine import QuerySet

from nolabs.infrastructure.settings import Settings
from nolabs.refined.domain.event_dispatcher import EventDispatcher
from nolabs.refined.domain.models import LocalisationJob, AminoAcid, Protein
from nolabs.refined.domain.models.experiment import Experiment


class Repository:
    _settings: Settings
    _event_dispatcher: EventDispatcher

    def __init__(self, settings: Settings, event_dispatcher: EventDispatcher):
        if not settings:
            raise ValueError("Settings")

        self._settings = settings
        self._event_dispatcher = event_dispatcher

    @property
    def experiments(self) -> QuerySet:
        return Experiment.objects

    @property
    def localisation_jobs(self) -> QuerySet:
        return LocalisationJob.objects

    @property
    def amino_acids(self) -> QuerySet:
        return AminoAcid.objects

    @property
    def proteins(self) -> QuerySet:
        return Protein.objects

