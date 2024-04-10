import { defineStore } from 'pinia';
import {ChemBLData, RcsbPdbData} from "src/api/client";

interface ChemblFunctionCallback  {
  (data: {files: ChemBLData[]}): void;
}

interface RcsbFunctionCallback {
  (data: {files: RcsbPdbData[]}): void;
}

interface StorageState {
  queryChemblEventHandlers: ChemblFunctionCallback[];
  queryRcsbPdbEventHandlers: RcsbFunctionCallback[];
}

export const useBioBuddyStore = defineStore('storage', {
  state: (): StorageState => ({
    queryChemblEventHandlers: [],
    queryRcsbPdbEventHandlers: [],
  }),
  actions: {
    addQueryChemblEventHandler(callback: ChemblFunctionCallback) {
      this.queryChemblEventHandlers.push(callback);
    },
    addQueryRcsbPdbEventHandler(callback: RcsbFunctionCallback) {
      this.queryRcsbPdbEventHandlers.push(callback);
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
  },
});
