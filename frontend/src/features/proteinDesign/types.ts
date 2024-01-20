export type ExperimentListItem = {
  id: string;
  name: string;
} | null;

export type Experiment = {
  id: string;
  name: string;
  pdbsContent: string[];
} | null;

export type InferenceRequest = {
  pdbFile: Blob;
  contig: string;
  numberOfDesigns: number;
  timesteps: number;
  hotspots: string;
  experimentName: string;
  experimentId?: string;
}