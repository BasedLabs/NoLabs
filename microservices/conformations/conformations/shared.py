import json
from enum import Enum
from multiprocessing import shared_memory, Lock

class SharedMemoryManager:
    def __init__(self, name='job_dict', size=1024):
        try:
            self.shm = shared_memory.SharedMemory(name=name, create=True, size=size)
            self.shm.buf[:size] = bytearray(size)  # Initialize shared memory with zeros
        except FileExistsError:
            self.shm = shared_memory.SharedMemory(name=name)
        self.size = size
        self.lock = Lock()
        self._initialize_dict()

    def _initialize_dict(self):
        if self.shm.buf[:4].tobytes() == b'\x00\x00\x00\x00':
            self._write_dict({})

    def _read_dict(self):
        with self.lock:
            raw_data = self.shm.buf[:self.size].tobytes().rstrip(b'\x00')
            if raw_data:
                return json.loads(raw_data.decode())
            return {}

    def _write_dict(self, data):
        with self.lock:
            serialized_data = json.dumps(data).encode()
            self.shm.buf[:len(serialized_data)] = serialized_data
            self.shm.buf[len(serialized_data):] = b'\x00' * (self.size - len(serialized_data))

    def check_job_exists(self, job_id: str) -> bool:
        data = self._read_dict()
        return job_id in data

    def get_job_status(self, job_id: str) -> 'JobStatusEnum':
        data = self._read_dict()
        return data.get(job_id) or JobStatusEnum.idle

    def change_status(self, job_id: str, status: 'JobStatusEnum'):
        data = self._read_dict()
        data[job_id] = status
        self._write_dict(data)

    def __del__(self):
        self.shm.close()
        self.shm.unlink()


class JobStatusEnum(str, Enum):
    running = 'running'
    idle = 'idle'


memory_manager = SharedMemoryManager()
