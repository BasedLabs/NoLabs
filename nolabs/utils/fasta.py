from typing import Dict, List

from nolabs.domain.amino_acid import AminoAcid


class FastaReader:
    def __init__(self):
        pass

    def get_ids2seqs(self, fasta_contents: str | List[str]) -> List[AminoAcid]:
        """
        Returns a dictionary {'sequence_id' : 'sequence'}
        """
        sequence_ids = []
        sequences = []
        current_sequence_id = None
        sequence_data = []

        if isinstance(fasta_contents, str):
            fasta_contents = [l for l in [l.strip() for l in fasta_contents.split('\n')] if l]

        for line in fasta_contents:
            if line.startswith('>'):
                if current_sequence_id:
                    sequences.append(''.join(sequence_data))
                    sequence_data = []

                current_sequence_id = line[1:].strip()
                sequence_ids.append(current_sequence_id)
            else:
                sequence_data.append(line.strip())

        if current_sequence_id:
            sequences.append(''.join(sequence_data))

        return [AminoAcid(name=sequence_id.split('|')[0], sequence=sequence) for sequence_id, sequence in zip(sequence_ids, sequences)]

    def get_ids2seqs_from_path(self, path: str) -> Dict[str, str]:
        """
        Returns a dictionary {'sequence_id' : 'sequence'}
        """
        with open(path, 'r') as f:
            contents = f.read()
            return self.get_ids2seqs(contents)

    def get_contents_from_path(self, path: str) -> str:
        """
        Returns plain content of a fasta file
        """
        with open(path, 'r') as f:
            contents = f.read()
            return contents

class FastaWriter:
    def __init__(self):
        pass

    def write_single_fasta(self, sequence_id: str, sequence: str, file_path: str) -> None:
        with open(file_path, 'w') as f:
            contents = ">" + sequence_id + "\n" + sequence
            f.write(contents)

    def write_multiple_fasta(self, sequence_ids: List[str], sequences: List[str], file_path: str) -> None:
        with open(file_path, 'w') as f:
            for sequence_id, sequence in zip(sequence_ids, sequences):
                self.write_single_fasta(sequence_id, sequence, file_path)

    def save_to_file(self, fasta_contents: str, file_path: str):
        with open(file_path, 'w') as f:
            f.write(fasta_contents)


