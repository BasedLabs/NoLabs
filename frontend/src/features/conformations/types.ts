import type {IntegratorsRequest} from "src/api/client";
import {Timeline} from "src/components/types";

export type ExperimentProperties = {
  pdbFile: File | null,
  totalFrames: number,
  temperatureK: number,
  takeFrameEvery: number,
  stepSize: number,
  replaceNonStandardResidues: boolean,
  addMissingAtoms: boolean,
  addMissingHydrogens: boolean,
  frictionCoeff: number,
  ignoreMissingAtoms: boolean,
  integrator: IntegratorsRequest,
};



export type Experiment = {
  id: string;
  name: string;
  pdbContent: File | null;
  timeline: Timeline[],
  properties: ExperimentProperties
} | null;

export type InferenceRequest = {
  pdbFile: File,
  experimentName: string,
  experimentId: string,
  totalFrames: number,
  temperatureK: number,
  takeFrameEvery: number,
  stepSize: number,
  replaceNonStandardResidues: boolean,
  addMissingAtoms: boolean,
  addMissingHydrogens: boolean,
  frictionCoeff: number,
  ignoreMissingAtoms: boolean,
  integrator: IntegratorsRequest,
}