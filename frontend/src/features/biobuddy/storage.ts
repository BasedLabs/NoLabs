import { defineStore } from 'pinia';
import type {FunctionCallReturnData} from "../../api/client";

interface FunctionCallback {
  (data: FunctionCallReturnData): void;
}

interface StorageState {
  queryChemblEventHandlers: FunctionCallback[];
  queryRcsbPdbEventHandlers: FunctionCallback[];
}

export const useBioBuddyStore = defineStore('storage', {
  state: (): StorageState => ({
    queryChemblEventHandlers: [],
    queryRcsbPdbEventHandlers: [],
  }),
  actions: {
    addQueryChemblEventHandler(callback: FunctionCallback) {
      this.queryChemblEventHandlers.push(callback);
    },
    addQueryRcsbPdbEventHandler(callback: FunctionCallback) {
      this.queryRcsbPdbEventHandlers.push(callback);
    },
    invokeQueryChemblEventHandlers(data: any) {
      this.queryChemblEventHandlers.forEach(handler => handler(data));
    },
    invokeQueryRcsbPdbEventHandlers(data: any) {
      this.queryRcsbPdbEventHandlers.forEach(handler => handler(data));
    },
  },
});
