export type ExperimentListItem = {
  id: string;
  name: string;
} | null;

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

export const PdbViews = {
  default: { key: 'default', title: 'Default representation' },
  cartoon: { key: 'cartoon', title: 'Cartoon' },
  backbone: { key: 'backbone', title: 'Backbone' },
  ballsAndSticks: { key: 'ball+stick', title: 'Balls and sticks' },
  contact: { key: 'contact', title: 'Contact' },
  helixorient: { key: 'helixorient', title: 'Helixorient' },
  hyperball: { key: 'hyperball', title: 'Hyperball' },
  licorice: { key: 'licorice', title: 'Licorice' },
  ribbon: { key: 'ribbon', title: 'Ribbon' },
  rope: { key: 'rope', title: 'Rope' },
  surface: { key: 'surface', title: 'Surface' },
  spacefill: { key: 'spacefill', title: 'Spacefill' },
  unitcell: { key: 'unitcell', title: 'Unitcell' }
};
