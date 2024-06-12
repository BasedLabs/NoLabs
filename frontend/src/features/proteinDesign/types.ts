export type JobProperties = {
  inputPdbFile: File | null;
  contig: string;
  numberOfDesigns: number;
  timesteps: number;
  hotspots: string;
};

export type Job = {
  id: string;
  name: string;
  generatedPdbs: File[];
  properties: JobProperties
} | null;

export type InferenceRequest = {
  pdbFile: File;
  contig: string;
  numberOfDesigns: number;
  timesteps: number;
  hotspots: string;
  jobName: string;
  jobId?: string;
}

