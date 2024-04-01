interface Smiles {
  smiles: string;
  drugLikeness: number;
  score: number;
  createdAt: Date;
}

interface Params {
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

interface Job {
  id: string;
  name: string;
  createdAt: string;
  running: boolean;
  learningCompleted: boolean;
}

interface Logs {
  output: string;
  dockingOutput: string;
  errors: string;
}
