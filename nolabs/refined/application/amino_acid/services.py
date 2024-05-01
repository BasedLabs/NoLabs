from typing import List

from nolabs.refined.application.amino_acid.api_models import RunAminoAcidRequest
from nolabs.refined.domain.models import Experiment, Protein, ProteinName
from nolabs.utils.fasta import FastaReader


async def get_input_proteins(experiment: Experiment, request: RunAminoAcidRequest) -> List[Protein]:
    input_amino_acids: List[Protein] = []
    fasta_reader = FastaReader()
    for fasta in request.fastas:
        fasta_content = await fasta.read()
        for amino_acid in fasta_reader.get_ids2seqs(fasta_content.decode('utf-8')):
            amino_acid = Protein.create(experiment=experiment, name=ProteinName(amino_acid.name),
                                        fasta_content=amino_acid.sequence)
            input_amino_acids.append(amino_acid)
            amino_acid.save()
    return input_amino_acids
