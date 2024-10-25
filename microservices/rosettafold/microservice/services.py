import glob
import math
import multiprocessing
import os
import pathlib
import subprocess
import threading
import uuid
from contextlib import contextmanager

import psutil

from microservice.api_models import RunRosettaFoldResponse
from microservice.loggers import logger

from leaf import ObjectType, FileObject, DirectoryObject


class RosettaService:
    _store = set()

    def __init__(self):
        self._root = DirectoryObject('/RoseTTAFold')
        self._fasta_name = '_run_' + str(uuid.uuid4())
        self._fasta_file_name = self._fasta_name + '.fasta'

    @contextmanager
    @staticmethod
    def _track_job(job_id: str | None):
        if job_id:
            RosettaService._store.add(job_id)
        try:
            yield True
        finally:
            if job_id:
                RosettaService._store.remove(job_id)

    def is_job_running(self, job_id: str) -> bool:
        return job_id in RosettaService._store

    async def run_rosettafold(self, job_id: str, fasta: bytes | None, a3m: bytes | None) -> RunRosettaFoldResponse:
        try:
            with RosettaService._track_job(job_id):
                return await self._run_rosettafold(fasta, a3m)
        except:
            logger.rosetta_exception()
            return RunRosettaFoldResponse(
                errors=['Critical error. Open the issue on our github and attach the error log'],
                pdb_content=None)
        finally:
            self.cleanup_generated_files()

    def cleanup_generated_files(self):
        for directory in self._root.directories.where(lambda o: '_run_' in o.name):
            directory.delete()

    async def _run_rosettafold(self, fasta: bytes, a3m: bytes) -> RunRosettaFoldResponse:
        bfd = self._root.directories.first_or_default(lambda o: o.name == 'bfd')
        uniref = self._root.directories.first_or_default(lambda o: o.name == 'UniRef30_2020_06')
        pdb100 = self._root.directories.first_or_default(lambda o: o.name == 'pdb100_2021Mar03')

        errors = []
        if not bfd:
            logger.bfd_directory_does_not_exist()
            errors.append('BFD directory does not exist. Check readme')

        if not uniref:
            logger.uniref_directory_does_not_exist()
            errors.append('Uniref directory does not exist. Check readme')

        if not pdb100:
            logger.pdb100_directory_does_not_exist()
            errors.append('PDB100 directory does not exist. Check readme')

        if errors:
            return RunRosettaFoldResponse(pdb_content=None, errors=errors)

        cpu_count = multiprocessing.cpu_count()
        gb = math.floor(psutil.virtual_memory().total / (1 << 30))

        if a3m:
            models_dir = self._root.add_directory(self._fasta_name).add_directory('models')
            msa_file = self._root.add_directory(self._fasta_name).add_file(self._fasta_name + 'msa0.a3m')
            msa_file.write_bytes(a3m)
            script = f'./run_RF2_msa.sh {msa_file.full_path} {os.path.join(models_dir.full_path, 'model')}'
        else:
            sequence = fasta.decode()
            if '>' not in sequence:
                sequence = f'>sequence\n{sequence}'
            self._root.add_file(self._fasta_file_name).write_bytes(sequence.encode())
            runner = self._root.files.first_or_default(lambda o: o.name == 'run_RF2.sh')
            script = f'CPU={cpu_count} MEM={gb} ./{runner.name} {self._fasta_file_name} -o {self._fasta_name}'

        logger.log_arbitrary(f'Running: {script}')

        def read_stdout(pipe):
            for line in iter(pipe.readline, ''):
                logger.rosetta_stdout(line)
            pipe.close()

        def read_stderr(pipe):
            for line in iter(pipe.readline, ''):
                logger.rosetta_stderr(line)
            pipe.close()

        process = subprocess.Popen(script, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True,
                                   cwd=self._root.full_path)

        stdout_thread = threading.Thread(target=read_stdout, args=(process.stdout,))
        stderr_thread = threading.Thread(target=read_stderr, args=(process.stderr,))

        stdout_thread.start()
        stderr_thread.start()

        stdout_thread.join()
        stderr_thread.join()

        return_code = process.wait()
        logger.rosetta_return_code(return_code)

        working_dir = self._root.directories.first_or_default(lambda o: o.name == self._fasta_name)
        if working_dir:
            pdb = working_dir.files.first_or_default(lambda o: o.name == 'model_00_pred.pdb', recursive=True)
            if pdb:
                return RunRosettaFoldResponse(pdb_content=pdb.read_string(), errors=[])

        logger.rosetta_fatal_error()
        return RunRosettaFoldResponse(pdb_content=None,
                                      errors=['Critical error. Open the issue on our github and attach the error log'])
