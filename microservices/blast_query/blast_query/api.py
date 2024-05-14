import os

from fastapi import FastAPI, HTTPException
from blast_query.api_models import SequenceQuery
from Bio.Blast import NCBIWWW
import xmltodict

app = FastAPI(
    title="BLAST Query"
)

def choose_dataset(program):
    if program in ["blastn", "tblastx", 'tblastn']:
        return "nt"
    elif program in ["blastp", "blastx"]:
        return "nr"
    else:
        raise HTTPException(status_code=400, detail="Invalid program")

def run_blast(program, query, descriptions=10, alignments=10, hitlist_size=10, expect=10.0):
    NCBIWWW.email = os.getenv("EMAIL")

    try:
        dataset = choose_dataset(program)
        result_handle = NCBIWWW.qblast(program, dataset, query, descriptions=descriptions, alignments=alignments, hitlist_size=hitlist_size, expect=expect)
        result_dict = xmltodict.parse(result_handle.read())
        result_handle.close()
        return result_dict
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/blastn")
async def blastn(query: SequenceQuery):
    """BLASTN compares a nucleotide query sequence against a nucleotide database.

    BLASTN is used for comparing a nucleotide sequence (DNA or RNA) against a nucleotide database. This is useful for finding regions of local similarity, which can suggest functional and evolutionary relationships between the sequences. BLASTN can align entire sequences or just portions that are similar. It’s particularly useful for identifying species, finding genes in a newly sequenced DNA, or mapping sequences to genomes.
    """
    return run_blast("blastn", query.sequence, query.descriptions, query.alignments, query.hitlist_size, query.expect)

@app.post("/blastp")
async def blastp(query: SequenceQuery):
    """BLASTP compares an amino acid (protein) query sequence against a protein database.

    BLASTP is used to compare protein sequences to sequence databases and calculates the statistical significance of matches. BLASTP can identify probable function and evolutionary relationships of a protein by locating matches to known protein sequences. Because proteins are directly responsible for function in living organisms, identifying and understanding protein sequences is crucial for biological research.
    """
    return run_blast("blastp", query.sequence, query.descriptions, query.alignments, query.hitlist_size, query.expect)

@app.post("/blastx")
async def blastx(query: SequenceQuery):
    """BLASTX compares a nucleotide sequence translated in all reading frames to a protein sequence database.

    BLASTX takes a nucleotide sequence, translates it in all six possible reading frames (three frames in each direction), and compares it against a protein database. This is useful for identifying potential protein products of an uncharacterized nucleotide sequence, checking for the presence of a protein in different species, or identifying potential protein-coding regions in a DNA sequence.
    """
    return run_blast("blastx", query.sequence, query.descriptions, query.alignments, query.hitlist_size, query.expect)

@app.post("/tblastn")
async def tblastn(query: SequenceQuery):
    """TBLASTN compares a protein sequence against a nucleotide sequence database that is dynamically translated into all reading frames.

    TBLASTN is used to find protein sequences that have been encoded in the nucleotide sequences present in a database. This tool translates the nucleotide database into all six possible reading frames for comparison with a protein query. This is particularly useful for identifying protein sequences in a genome where the gene might not have been previously identified or annotated.
    """
    return run_blast("tblastn", query.sequence, query.descriptions, query.alignments, query.hitlist_size, query.expect)

@app.post("/tblastx")
async def tblastx(query: SequenceQuery):
    """TBLASTX compares the six-frame translations of a nucleotide query against the six-frame translations of a nucleotide sequence database.

    TBLASTX is a computationally intensive tool that translates a nucleotide query sequence in all six frames and compares it against a nucleotide database that is also translated in all six frames. This can be useful for comparing different genomic regions or whole genomes to find potential coding regions across species that do not have annotated proteins. It provides a way to identify conserved protein sequences in organisms that are too divergent for nucleotide-level comparisons to be effective.
    """
    return run_blast("tblastx", query.sequence, query.descriptions, query.alignments, query.hitlist_size, query.expect)
