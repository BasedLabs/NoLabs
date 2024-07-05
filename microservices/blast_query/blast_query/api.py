import os
from typing import Any, Dict

from fastapi import FastAPI, HTTPException
from blast_query.api_models import SequenceQuery, BlastType
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
    raise HTTPException(status_code=400, detail="Invalid program")

def run_blast(program: BlastType, query: str, descriptions: int = 10, alignments: int = 10, hitlist_size: int = 10, expect: float = 10.0) -> Dict[str, Any]:
    NCBIWWW.email = os.getenv("EMAIL")

    try:
        dataset = choose_dataset(program)
        result_handle = NCBIWWW.qblast(program.value, dataset, query, descriptions=descriptions, alignments=alignments, hitlist_size=hitlist_size, expect=expect)
        result_dict = xmltodict.parse(result_handle.read())
        result_handle.close()
        return result_dict
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/blast")
async def blast(query: SequenceQuery) -> Dict[str, Any]:
    """
    Perform a BLAST query with the specified type.

    Args:
        query (SequenceQuery): The query parameters including the sequence and BLAST type.

    Returns:
        Dict[str, Any]: The result of the BLAST query.
    """
    return run_blast(query.type, query.sequence, query.descriptions, query.alignments, query.hitlist_size, query.expect)
