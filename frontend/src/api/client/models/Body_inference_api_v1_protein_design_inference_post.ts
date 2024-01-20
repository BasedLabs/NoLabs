/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
export type Body_inference_api_v1_protein_design_inference_post = {
    experiment_name: string;
    experiment_id?: string;
    pdb_file: Blob;
    contig?: string;
    number_of_desings?: number;
    timesteps?: number;
    hotspots?: string;
};

