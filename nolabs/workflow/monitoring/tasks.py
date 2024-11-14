from typing import Iterable, Tuple, Dict, List

from nolabs.infrastructure.redis_client_factory import rd


class OrphanedTasksTracker:
    _keys_pattern = "orphaned_tasks_tracker:"

    @classmethod
    def track(cls, task_id: str, queue: str, timestamp: float):
        rd.hset(f"{cls._keys_pattern}{task_id}", mapping={
            'task_id': task_id,
            'queue': queue,
            'timestamp': timestamp
        })

    @classmethod
    def get_task_data(cls) -> List[Dict[str, str]]:
        result = []
        for key in rd.keys(f"{cls._keys_pattern}*"):
            task_data = rd.hgetall(key)
            if task_data:
                task_data['timestamp'] = float(task_data['timestamp'])
                result.append(task_data)
        return result

    @classmethod
    def remove_task(cls, task_id):
        key = f"{cls._keys_pattern}{task_id}"
        rd.delete(key)
