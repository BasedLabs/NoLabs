import { defineStore } from 'pinia';
import { ChemBLData, RcsbPdbData } from "src/refinedApi/client";

interface ChemblFunctionCallback {
  (data: { files: ChemBLData[] }): void;
}

interface RcsbFunctionCallback {
  (data: { files: RcsbPdbData[] }): void;
}

interface WorkflowAdjustmentCallback {
  (workflowData: any): void;
}

interface StorageState {
  queryChemblEventHandlers: ChemblFunctionCallback[];
  queryRcsbPdbEventHandlers: RcsbFunctionCallback[];
  workflowAdjustmentEventHandlers: WorkflowAdjustmentCallback[];
}

export const useBioBuddyStore = defineStore('storage', {
  state: (): StorageState => ({
    queryChemblEventHandlers: [],
    queryRcsbPdbEventHandlers: [],
    workflowAdjustmentEventHandlers: [],
  }),
  actions: {
    addQueryChemblEventHandler(callback: ChemblFunctionCallback) {
      this.queryChemblEventHandlers.push(callback);
    },
    addQueryRcsbPdbEventHandler(callback: RcsbFunctionCallback) {
      this.queryRcsbPdbEventHandlers.push(callback);
    },
    addWorkflowAdjustmentEventHandler(callback: WorkflowAdjustmentCallback) {
      this.workflowAdjustmentEventHandlers.push(callback);
    },
    async invokeQueryChemblEventHandlers(data: any) {
      for (const handler of this.queryChemblEventHandlers) {
        await handler(data);
      }
    },
    async invokeQueryRcsbPdbEventHandlers(data: any) {
      for (const handler of this.queryRcsbPdbEventHandlers) {
        await handler(data);
      }
    },
    async invokeWorkflowAdjustment(data: any) {
      for (const handler of this.workflowAdjustmentEventHandlers) {
        await handler(data);
      }
    },
  },
});
