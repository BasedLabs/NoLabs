from typing import List

from nolabs.refined.application.controllers.amino_acid.api_models import RunAminoAcidRequest
from nolabs.refined.domain.models import Experiment, AminoAcid, AminoAcidName
from nolabs.utils.fasta import FastaReader


async def get_input_amino_acids(experiment: Experiment, request: RunAminoAcidRequest) -> List[AminoAcid]:
    input_amino_acids: List[AminoAcid] = []
    fasta_reader = FastaReader()
    for fasta in request.fastas:
        fasta_content = await fasta.read()
        for amino_acid in fasta_reader.get_ids2seqs(fasta_content.decode('utf-8')):
            amino_acid = AminoAcid.create(experiment=experiment, name=AminoAcidName(amino_acid.name),
                                          content=amino_acid.sequence)
            input_amino_acids.append(amino_acid)
            amino_acid.save()
    return input_amino_acids