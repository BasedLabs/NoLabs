export type ExperimentProperties = {
  aminoAcidSequence: string | null | undefined,
  fastas: Array<File>
};

export type AminoAcid = {
  name: string;
  sequence: string;
  cytosolicProteins: number;
  mitochondialProteins: number;
  nuclearProteins: number;
  otherProteins: number;
  extracellularSecretedProteins: number;
};

export type Experiment = {
  id: string;
  name: string;
  aminoAcids: Array<AminoAcid>,
  properties: ExperimentProperties
} | null;

export type InferenceRequest = {
  experimentName: string;
  experimentId?: string;
  aminoAcidSequence: string | undefined,
  fastas: Array<File>
}

