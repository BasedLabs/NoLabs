/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ProteinLocalisationResponse } from './ProteinLocalisationResponse';
export type ProteinResponse = {
    id: string;
    name: string;
    experiment_id: string;
    binding_pockets: Array<number>;
    fasta_name: string;
    pdb_name: string;
    localisation?: (ProteinLocalisationResponse | null);
    gene_ontology?: (Record<string, any> | null);
    soluble_probability?: (number | null);
    msa?: (string | null);
    md_pdb_content?: (string | null);
    fasta_content?: (string | null);
    pdb_content?: (string | null);
};

