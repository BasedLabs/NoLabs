import io
from typing import List
from Bio import SeqIO
from werkzeug.datastructures import FileStorage


def get_sequences(files: List[FileStorage]):
    for file in files:
        content = file.read().decode('utf-8')
        fasta_sequences = SeqIO.parse(io.StringIO(content), 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            yield sequence
