from typing import List, Iterable

from nolabs.infrastructure.redis_client_factory import Redis


class OrphanedTasksTracker:
    def __init__(self):
        self._keys_pattern = 'orphaned_tasks_tracker:'

    def track(self, task_id):
        Redis.client.set(f'{self._keys_pattern}{task_id}', 1)

    def get_task_ids(self) -> Iterable[str]:
        for key in Redis.client.keys(f'{self._keys_pattern}*'):
            yield key.replace(self._keys_pattern, '')

    def remove_task(self, task_id):
        key = f'{self._keys_pattern}{task_id}'
        Redis.client.delete(key)
