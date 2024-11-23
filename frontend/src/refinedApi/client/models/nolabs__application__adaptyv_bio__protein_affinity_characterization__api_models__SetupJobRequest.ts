/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
export type nolabs__application__adaptyv_bio__protein_affinity_characterization__api_models__SetupJobRequest = {
    job_id: string;
    /**
     * Every protein counts as one design. The more designs you test the cheaper it becomes.
     */
    number_of_designs: number;
    /**
     * AA (avg. protein length)
     */
    dna_length: number;
    /**
     * A replicate is a single test for a protein design. The more replicates you select the more confidence you can have in your data.
     */
    replicates: number;
    /**
     * Email address to send the report
     */
    report_email: string;
    /**
     * Target id of the experiment
     */
    target_id: string;
    /**
     * Total price of the experiment
     */
    cart_total: number;
    /**
     * URL of the session where the experiment is located
     */
    session_url: string;
    swissprot_id: string;
};

