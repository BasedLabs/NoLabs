# FastAPI BLAST Service

This project provides a simple FastAPI service that allows users to perform BLAST searches using various BLAST databases. It supports nucleotide and protein BLAST queries through different endpoints.

NCBI guidelines, from https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo state:
1. Do not contact the server more often than once every 10 seconds.
2. Do not poll for any single RID more often than once a minute.
3. Use the URL parameter email and tool, so that the NCBI can contact you if there is a problem.
4. Run scripts weekends or between 9 pm and 5 am Eastern time on weekdays if more than 50 searches
will be submitted.

Tool description:
1.  BLASTN compares a nucleotide query sequence against a nucleotide database.

    BLASTN is used for comparing a nucleotide sequence (DNA or RNA) against a nucleotide database. This is useful for finding regions of local similarity, which can suggest functional and evolutionary relationships between the sequences. BLASTN can align entire sequences or just portions that are similar. It’s particularly useful for identifying species, finding genes in a newly sequenced DNA, or mapping sequences to genomes.
2.  BLASTP compares an amino acid (protein) query sequence against a protein database.

    BLASTP is used to compare protein sequences to sequence databases and calculates the statistical significance of matches. BLASTP can identify probable function and evolutionary relationships of a protein by locating matches to known protein sequences. Because proteins are directly responsible for function in living organisms, identifying and understanding protein sequences is crucial for biological research.
3.  BLASTX compares a nucleotide sequence translated in all reading frames to a protein sequence database.

    BLASTX takes a nucleotide sequence, translates it in all six possible reading frames (three frames in each direction), and compares it against a protein database. This is useful for identifying potential protein products of an uncharacterized nucleotide sequence, checking for the presence of a protein in different species, or identifying potential protein-coding regions in a DNA sequence.
4.  TBLASTN compares a protein sequence against a nucleotide sequence database that is dynamically translated into all reading frames.

    TBLASTN is used to find protein sequences that have been encoded in the nucleotide sequences present in a database. This tool translates the nucleotide database into all six possible reading frames for comparison with a protein query. This is particularly useful for identifying protein sequences in a genome where the gene might not have been previously identified or annotated.
5.  TBLASTX compares the six-frame translations of a nucleotide query against the six-frame translations of a nucleotide sequence database.

    TBLASTX is a computationally intensive tool that translates a nucleotide query sequence in all six frames and compares it against a nucleotide database that is also translated in all six frames. This can be useful for comparing different genomic regions or whole genomes to find potential coding regions across species that do not have annotated proteins. It provides a way to identify conserved protein sequences in organisms that are too divergent for nucleotide-level comparisons to be effective.

Useful resources:
- BLAST documentation https://blast.ncbi.nlm.nih.gov/Blast.cgi
- Biopython BLAST documentation here: https://biopython.org/docs/1.75/api/Bio.Blast.html
- Biopython tutorial: https://biopython.org/DIST/docs/tutorial/Tutorial.pdf
