/* generated using openapi-typescript-codegen -- do no edit */
/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
export type nolabs__application__adaptyv_bio__protein_affinity_characterization__api_models__SetupJobRequest = {
    job_id: string;
    /**
     * Every protein counts as one design. The more designs you test the cheaper it becomes.
     */
    number_of_designs?: (number | null);
    /**
     * AA (avg. protein length)
     */
    dna_length?: (number | null);
    /**
     * A replicate is a single test for a protein design. The more replicates you select the more confidence you can have in your data.
     */
    replicates?: (number | null);
    /**
     * Email address to send the report
     */
    report_email?: (string | null);
    /**
     * Target id of the experiment
     */
    target_id?: (string | null);
    /**
     * Target name of the experiment
     */
    target_name?: (string | null);
    /**
     * Total price of the experiment
     */
    cart_total?: (number | null);
    swissprot_id?: (string | null);
};

