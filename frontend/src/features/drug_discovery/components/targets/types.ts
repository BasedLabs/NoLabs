import {TargetMetaData} from "src/api/client/models/TargetMetaData.ts";
import {LigandMetaData} from "src/api/client/models/LigandMetaData.ts";

export interface ExtendedTargetMetaData extends TargetMetaData {
  loadingLigands: boolean;
  ligands: ExtendedLigandMetaData[];
  loadingTargetData: boolean;
  data?: {
    proteinSequence?: string;
    pdbContents?: string | null;
  };
}

export interface ExtendedLigandMetaData extends LigandMetaData {
  loadingLigandData: boolean;
  jobs: any[]; // Define a more specific type if possible
  data?: {
    sdf_file?: string;
    smiles?: string;
  };
}
