export type AminoAcidJobProperties = {
    fastas: Array<File>
};

export type Job<TAminoAcid> = {
    id: string;
    name: string;
    aminoAcids: Array<TAminoAcid>,
    properties: AminoAcidJobProperties
} | null;

export type InferenceRequest = {
    jobName: string;
    jobId?: string;
    fastas: Array<File>
}
