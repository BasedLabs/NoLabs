/* generated using openapi-typescript-codegen -- do not edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ProteinLocalisationResponse } from './ProteinLocalisationResponse';
export type ProteinMetadataResponse = {
    id: string;
    name: string;
    experiment_id: string;
    binding_pockets: Array<number>;
    fasta_name: string;
    pdb_name: string;
    localisation?: (ProteinLocalisationResponse | null);
    gene_ontology?: (Record<string, any> | null);
    soluble_probability?: (number | null);
    link?: (string | null);
};

