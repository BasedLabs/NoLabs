interface Smiles {
  smiles: string;
  drugLikeness: number;
  score: number;
  stage: string;
  createdAt: Date;
}

interface Status {
  running: boolean;
  samplingAllowed: boolean;
}

interface Properties {
  pdbFile: File | null;
  centerX: number;
  centerY: number;
  centerZ: number;
  sizeX: number;
  sizeY: number;
  sizeZ: number;
  batchSize: number;
  minscore: number;
  epochs: number;
}

interface Experiment {
  id: string;
  name: string;
  running: boolean;
  samplingAllowed: boolean;
  properties: Properties;
  createdAt: Date;
}

interface Logs {
  output: string;
  dockingOutput: string;
  errors: string;
}
