import { defineStore } from 'pinia';
import { LigandResponse, ProteinResponse } from 'src/refinedApi/client';
import { deleteProtein, getAllProteins, uploadProtein, updateProteinName } from 'src/features/drug_discovery/refinedApi';

export const useWorkflowStore = defineStore('workflowStore', {
    state: () => ({
        proteins: [] as ProteinResponse[],
        ligands: [] as LigandResponse[],
    }),
    actions: {
        async getAllProteins(experimentId: string) {
            const response = await getAllProteins(experimentId);
            this.proteins = response;
        },
        async uploadProteinToExperiment(experimentId: string, name?: string, fastaFile?: Blob, pdbFile?: Blob, metaData?: Record<string, string>) {
            try {
                const uploadedProtein = await uploadProtein(
                    experimentId,
                    name,
                    fastaFile,
                    pdbFile
                );
                this.proteins.push(uploadedProtein);
            } catch (error) {
                console.error('Error uploading protein:', error);
            }
        },
        async deleteProteinFromExperiment(proteinId: string) {
            try {
                await deleteProtein(proteinId);
                this.proteins = this.proteins.filter(protein => protein.id !== proteinId);
            } catch (error) {
                console.error('Error deleting protein:', error);
            }
        },
        async changeProteinName(proteinId: string, newName: string) {
            try {
                await updateProteinName(proteinId, newName);
                const protein = this.proteins.find(protein => protein.id === proteinId);
                if (protein) {
                  protein.name = newName;
                }
            } catch (error) {
                console.error('Error deleting protein:', error);
            }
        }
    }
});
