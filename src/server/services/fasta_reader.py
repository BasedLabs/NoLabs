import io

from Bio import SeqIO


def get_sequences(files):
    for file in files:
        content = file.read().decode('utf-8')
        fasta_sequences = SeqIO.parse(io.StringIO(content), 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            yield sequence
