import {AminoAcid} from "src/features/aminoAcid/localisation/types";

export type AminoAcidExperimentProperties = {
    aminoAcidSequence: string | null | undefined,
    fastas: Array<File>
};

export type Experiment<TAminoAcid> = {
    id: string;
    name: string;
    aminoAcids: Array<TAminoAcid>,
    properties: AminoAcidExperimentProperties
} | null;

export type InferenceRequest = {
    experimentName: string;
    experimentId?: string;
    aminoAcidSequence: string | undefined,
    fastas: Array<File>
}