import ssl
import urllib.request
from typing import Dict, Any

import xmltodict
from Bio.Blast import NCBIWWW

from api_models import BlastType


def choose_dataset(program):
    if program in ["blastn", "tblastx", 'tblastn']:
        return "nt"
    elif program in ["blastp", "blastx"]:
        return "nr"
    raise ValueError("Dataset is not defined")


def run_blast(program: BlastType,
              query: str,
              descriptions: int = 10,
              alignments: int = 10,
              hitlist_size: int = 10,
              expect: float = 10.0) -> Dict[str, Any]:
    ssl_context = ssl.create_default_context()
    ssl_context.check_hostname = False
    ssl_context.verify_mode = ssl.CERT_NONE
    dataset = choose_dataset(program)
    opener = urllib.request.build_opener(urllib.request.HTTPSHandler(context=ssl_context))
    urllib.request.install_opener(opener)
    result_handle = NCBIWWW.qblast(program.value, dataset, query, descriptions=descriptions, alignments=alignments,
                                   hitlist_size=hitlist_size, expect=expect)
    result_dict = xmltodict.parse(result_handle.read())
    result_handle.close()
    return result_dict
