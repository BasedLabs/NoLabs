class WsRequest:
    def __init__(self, amino_acid_seq, amino_acid_seq_files):
        self.amino_acid_seq = amino_acid_seq
        self.amino_acid_seq_files = amino_acid_seq_files

    @staticmethod
    def from_client_json(json):
        return WsRequest(json['aminoAcidSeq'], json['aminoAcidSeqFiles'])