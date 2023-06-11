import io

from Bio import SeqIO


def get_sequences(files):
    for file in files:
        fasta_sequences = SeqIO.parse(io.StringIO(file), 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            yield sequence
