export type ExperimentProperties = {
  inputPdbFile: File | null;
  contig: string;
  numberOfDesigns: number;
  timesteps: number;
  hotspots: string;
};

export type Experiment = {
  id: string;
  name: string;
  generatedPdbs: File[];
  properties: ExperimentProperties
} | null;

export type InferenceRequest = {
  pdbFile: File;
  contig: string;
  numberOfDesigns: number;
  timesteps: number;
  hotspots: string;
  experimentName: string;
  experimentId?: string;
}

