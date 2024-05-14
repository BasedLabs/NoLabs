# FastAPI BLAST Service

This project provides a simple FastAPI service that allows users to perform BLAST searches using various BLAST databases. It supports nucleotide and protein BLAST queries through different endpoints.

NCBI guidelines, from https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo state:
Do not contact the server more often than once every 10 seconds.
2. Do not poll for any single RID more often than once a minute.
3. Use the URL parameter email and tool, so that the NCBI can contact you if there is a problem.
4. Run scripts weekends or between 9 pm and 5 am Eastern time on weekdays if more than 50 searches
will be submitted.

Useful resources:
- BLAST documentation https://blast.ncbi.nlm.nih.gov/Blast.cgi
- Biopython BLAST documentation here: https://biopython.org/docs/1.75/api/Bio.Blast.html
- Biopython tutorial: https://biopython.org/DIST/docs/tutorial/Tutorial.pdf
