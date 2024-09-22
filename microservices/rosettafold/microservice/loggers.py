import logging
import sys

from fastapi import UploadFile
from microservice.api_models import RunRosettaFoldResponse
from pythonjsonlogger import jsonlogger

_logger = logging.getLogger()

_logger.setLevel(logging.DEBUG)

_logHandler = logging.StreamHandler(sys.stdout)
_formatter = jsonlogger.JsonFormatter()
_logHandler.setFormatter(_formatter)
_logger.addHandler(_logHandler)


class Log:
    def run_rosettafold_request(
        self, job_id: str, fasta: bytes | None, a3m: bytes | None
    ):
        _logger.info(
            "Run rosettafold request",
            extra={
                "job_id": job_id,
                "fasta": fasta.decode("utf-8")[:15] if fasta else None,
                "a3m": a3m.decode("utf-8")[:15] if a3m else None,
            },
        )

    def run_rosettafold_response(self, response: RunRosettaFoldResponse):
        d = response.as_log_dict()
        _logger.info("Run rosettafold response", extra=d)

    def bfd_directory_does_not_exist(self):
        _logger.error("BFD directory does not exist. Check readme")

    def uniref_directory_does_not_exist(self):
        _logger.error("Uniref directory does not exist. Check readme")

    def pdb100_directory_does_not_exist(self):
        _logger.error("PDB100 directory does not exist. Check readme")

    def rosetta_stdout(self, s: str):
        _logger.info("Rosetta stdout", extra={"stdout": s})

    def rosetta_stderr(self, s: str):
        _logger.error("Rosetta stderr", extra={"stderr": s})

    def rosetta_error(self, s: str, file: str):
        _logger.error("Rosetta error", extra={"content": s, "source_file": file})

    def rosetta_fatal_error(self):
        _logger.fatal(
            "Fatal error. Open the issue on our github and attach the error log"
        )

    def rosetta_return_code(self, rt: int):
        _logger.info("Rosetta finish", extra={"return_code": rt})

    def rosetta_exception(self):
        _logger.exception("Unexpected error. Open the issue on our github")

    def log_arbitrary(self, message: str):
        _logger.info(message)


logger = Log()
