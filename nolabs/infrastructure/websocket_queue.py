import json
from multiprocessing import shared_memory, Lock
from typing import Optional, Any, List, Dict


class WebsocketsQueue:
    def __init__(self, name='websockets_queue', size=1024 * 4):
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
            with self.lock:
                self._write_arr([])

    def read_last(self) -> Optional[Dict[str, Any]]:
        with self.lock:
            raw_data = self.shm.buf[:self.size].tobytes().rstrip(b'\x00').decode()
            if raw_data:
                queue = json.loads(raw_data)
                item = queue.pop()
                self._write_arr(queue)
                return item
            return None

    def write(self, data: Dict[str, Any]):
        with self.lock:
            queue: List[Dict[str, Any]] = json.loads(self.shm.buf[:self.size].tobytes().rstrip(b'\x00').decode())
            queue.append(data)
            self._write_arr(queue)

    def _write_arr(self, queue: List[Dict[str, Any]]):
        serialized_data = json.dumps(queue).encode()
        self.shm.buf[:len(serialized_data)] = serialized_data
        self.shm.buf[len(serialized_data):] = b'\x00' * (self.size - len(serialized_data))

    def __del__(self):
        self.shm.close()
        self.shm.unlink()


websockets_queue = WebsocketsQueue()
