from typing import Dict, List


class FastaReader:
    def __init__(self):
        pass

    def get_ids_and_sequences(self, fasta_contents: str) -> Dict[str, str]:
        """
        Returns a dictionary {'sequence_id' : 'sequence'}
        """
        sequence_ids = []
        sequences = []
        current_sequence_id = None
        sequence_data = []

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

        return {sequence_id: sequence for sequence_id, sequence in zip(sequence_ids, sequences)}

    def get_data_from_path(self, path: str) -> Dict[str, str]:
        """
        Returns a dictionary {'sequence_id' : 'sequence'}
        """
        with open(path, 'r') as f:
            contents = f.read()
            return self.get_ids_and_sequences(contents)

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


