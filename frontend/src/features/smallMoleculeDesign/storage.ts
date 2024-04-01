import {defineStore} from "pinia";


const useSmallMoleculesDesignStore = defineStore("smallMoleculesDesignStore", {
  actions: {
    async job(jobId: string): Promise<Job> {
      return {
        id: '1',
        name: 'Test 1',
        createdAt: new Date().toISOString(),
        running: false,
        learningCompleted: false
      };
    },
    async smilesData(jobId: string): Promise<Smiles[]>{
      return [];
    },
    async jobs(): Promise<Job[]> {
      return [
        {
          id: '1',
          name: 'Test 1',
          createdAt: new Date().toISOString(),
          running: false,
          learningCompleted: false
        },
        {
          id: '2',
          name: 'Test 2',
          createdAt: new Date().toISOString(),
          running: true,
          learningCompleted: false
        }
      ]
    },
    async logs(jobId: string): Promise<Logs> {
      return {
        output: 'asd1',
        dockingOutput: 'asd2',
        errors: 'asd3'
      };
    },
    async params(jobId: string): Promise<Params | null> {
      return null;
    },
    async startJob(jobId: string) {

    },
    async generateMoreMolecules(jobId: string) {

    },
    async stopJob(jobId: string) {

    },
    async deleteJob(jobId: string) {

    },
    async saveParams(jobId: string, params: Params): Promise<void> {

    }
  }
});

export default useSmallMoleculesDesignStore;
