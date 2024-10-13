from dotenv import load_dotenv

load_dotenv("infrastructure/.env")

from .mongo_connector import mongo_connect
from infrastructure.socket_server import start_queue_worker
mongo_connect()
start_queue_worker()