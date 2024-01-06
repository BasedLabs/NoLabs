from typing import Dict, List

import typing
from conformations_microservice import ForceFields, WaterForceFields
from fastapi import WebSocket


def _create_ws_message(message: str = None, error: str = None) -> Dict[str, str]:
    return {
        'message': message,
        'error': error
    }


class ConformationsWebsocket:
    def __init__(self, websocket: WebSocket | None):
        self._websocket = websocket

    def pdb_fixer_start(self):
        self._websocket.send_json(data=_create_ws_message('Pdb fixer started'))

    def pdb_fixer_finish(self):
        self._websocket.send_json(data=_create_ws_message())

    def pdb_fixer_errors(self, errors: List[str]):
        self._websocket.send_json(data=_create_ws_message(message='Pdb fixer errors', error=', '.join(errors)))

    def gromacs_simulations_start(self, force_field: ForceFields, water_force_fields: WaterForceFields):
        self._websocket.send_json(data=_create_ws_message(
            f'Gromacs simulations started, force field: {force_field}, water force field: {water_force_fields}'))

    def gromacs_simulations_finish(self):
        self._websocket.send_json(data=_create_ws_message('Gromacs simulations finished'))

    def gromacs_simulations_errors(self, errors: List[str]):
        self._websocket.send_json(
            data=_create_ws_message(message='Gromacs simulations errors', error=', '.join(errors)))

    def gromacs_gro_top_generator_start(self, force_field: ForceFields, water_force_fields: WaterForceFields):
        self._websocket.send_json(data=_create_ws_message(
            f'Gromacs gro top generator start, force field: {force_field}, water force field: {water_force_fields}'))

    def gromacs_gro_top_generator_finish(self):
        self._websocket.send_json(data=_create_ws_message('Gromacs gro top generator finish'))

    def gromacs_gro_top_generator_errors(self, errors: List[str]):
        self._websocket.send_json(
            data=_create_ws_message(message='Gromacs gro top generator errors', error=', '.join(errors)))

    def pdb_simulations_start(self, force_field: ForceFields, water_force_fields: WaterForceFields):
        self._websocket.send_json(data=_create_ws_message(
            f'Pdb simulations started, force field: {force_field}, water force field: {water_force_fields}'))

    def pdb_simulations_finish(self):
        self._websocket.send_json(data=_create_ws_message('Pdb simulations finished'))

    def pdb_simulations_errors(self, errors: List[str]):
        self._websocket.send_json(data=_create_ws_message(message='Pdb simulations errors', error=', '.join(errors)))
