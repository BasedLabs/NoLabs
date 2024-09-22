import asyncio
import functools
import ssl
import urllib.request
from enum import Enum
from typing import Any, Dict

import xmltodict
from Bio.Blast import NCBIWWW
from domain.exceptions import ErrorCodes, NoLabsException

from nolabs.infrastructure.settings import settings


class BlastType(str, Enum):
    BLASTN = "blastn"
    BLASTP = "blastp"
    BLASTX = "blastx"
    TBLASTN = "tblastn"
    TBLASTX = "tblastx"


async def _async_qblast(
    program, dataset, sequence, descriptions, alignments, hitlist_size, expect
):
    loop = asyncio.get_event_loop()

    func = functools.partial(
        NCBIWWW.qblast,
        program=program,
        database=dataset,
        sequence=sequence,
        descriptions=descriptions,
        alignments=alignments,
        hitlist_size=hitlist_size,
        expect=expect,
    )

    return await loop.run_in_executor(None, func)


async def run_blast(
    program: BlastType,
    sequence: str,
    descriptions: int = 10,
    alignments: int = 10,
    hitlist_size: int = 10,
    expect: float = 10.0,
) -> Dict[str, Any]:
    """
    :param sequence: Could be nucleotide sequence for blastn, tblastx, or tblastn, or amino acid sequence for blastp or blastx.
    """
    NCBIWWW.email = settings.blast_email

    # Create an SSL context that does not verify certificates
    ssl_context = ssl.create_default_context()
    ssl_context.check_hostname = False
    ssl_context.verify_mode = ssl.CERT_NONE

    try:
        datasets = {
            BlastType.BLASTN: "nt",
            BlastType.TBLASTX: "nt",
            BlastType.TBLASTN: "nt",
            BlastType.BLASTP: "nr",
            BlastType.BLASTX: "nr",
        }
        dataset = datasets[program]
        # Patch urllib to use the custom SSL context
        opener = urllib.request.build_opener(
            urllib.request.HTTPSHandler(context=ssl_context)
        )
        urllib.request.install_opener(opener)

        result_handle = await _async_qblast(
            program.value,
            dataset,
            sequence,
            descriptions=descriptions,
            alignments=alignments,
            hitlist_size=hitlist_size,
            expect=expect,
        )
        result_dict = xmltodict.parse(result_handle.read())
        result_handle.close()
        return result_dict
    except Exception as e:
        raise NoLabsException(ErrorCodes.blast_api_error) from e
