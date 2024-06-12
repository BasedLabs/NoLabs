import type {IntegratorsRequest} from "src/refinedApi/client";
import {Timeline} from "src/components/types";

export type JobProperties = {
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

export type Job = {
  id: string;
  name: string;
  pdbContent: File | null;
  timeline: Timeline[],
  properties: JobProperties
} | null;

export type InferenceRequest = {
  pdbFile: File,
  experimentId: string,
  jobName: string,
  jobId: string,
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
