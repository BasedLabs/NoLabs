<template>
    <q-list bordered separator class="bg-black">
        <q-item-label header class="text-white text-h6">Protein List</q-item-label>
        <q-item>
            <q-btn size="md" class="full-width q-pm-sm" push color="info" @click="uploadProteinDialog = true">
                + upload proteins
                <q-tooltip>Add proteins to the experiment</q-tooltip>
            </q-btn>
        </q-item>
        <q-expansion-item popup expand-separator :content-inset-level="1" v-for="protein in proteins" :key="protein.id"
            :label="protein.name">
            <q-btn color="negative" @click="deleteProtein(protein)" icon="delete" flat>
                Delete protein
            </q-btn>
            <ProteinDetail v-if="protein.fasta_content" :experiment-id="experimentId" :original-protein="protein">
            </ProteinDetail>
        </q-expansion-item>
    </q-list>
    <q-dialog v-model="uploadProteinDialog">
        <q-card>
            <q-card-section>
                <div class="text-h6">Upload Protein Files</div>
                <q-file v-model="proteinFileModel" ref="proteinFileInput" accept=".fasta" multiple
                    label="Choose files" />
            </q-card-section>

            <q-card-actions align="right">
                <q-btn flat label="Close" color="negative" @click="uploadProteinDialog = false" />
                <q-btn flat label="Upload" color="info" @click="uploadProteinFiles" />
            </q-card-actions>
        </q-card>
    </q-dialog>
</template>

<script lang="ts">
import { QSpinnerOrbit } from "quasar";
import { defineComponent } from "vue";
import { useWorkflowStore } from "src/features/drug_discovery/components/workflow/storage";
import { ProteinResponse } from "src/refinedApi/client";
import ProteinDetail from "./ProteinDetail.vue"
import { FunctionCallReturnData, RcsbPdbData } from "src/refinedApi/client";
import { useBioBuddyStore } from "src/features/biobuddy/storage";

export default defineComponent({
    name: "ProteinList",
    components: {
        ProteinDetail
    },
    props: {
        experimentId: {
            type: String,
            required: true,
        },
    },
    data() {
        return {
            uploadProteinDialog: false,
            proteinFileModel: [],
            selectedProtein: null,
            rscbQueryCallBack: null as ((data: FunctionCallReturnData) => void) | null,
            proteins: [] as ProteinResponse[]
        };
    },
    async mounted() {
        const workflowStore = useWorkflowStore();
        await workflowStore.getAllProteins(this.experimentId);
        this.proteins = workflowStore.proteins;


        const bioBuddyStore = useBioBuddyStore();
        this.rscbQueryCallBack = async (data: { files: RcsbPdbData[] }) => {
            const files: File[] = [];
            const metaDatasArray: Record<string, string>[] = [];

            for (let rcsbObject of data.files) {
                const structureId = this.extractStructureId(rcsbObject.metadata.link);
                const file = new File([new Blob([rcsbObject.content!])], structureId + ".fasta");
                files.push(file);
                const metadataAsString: Record<string, string> = {};
                Object.entries(rcsbObject.metadata).forEach(([key, value]) => {
                    metadataAsString[key] = typeof value === 'object' ? JSON.stringify(value) : String(value);
                });
                metaDatasArray.push(metadataAsString);
            }
            await this.handleProteinUpload(files, metaDatasArray);
        };
        bioBuddyStore.addQueryRcsbPdbEventHandler(this.rscbQueryCallBack);
    },
    methods: {
        uploadProteinFiles() {
            if (this.proteinFileModel.length > 0) {
                this.handleProteinUpload(this.proteinFileModel)
                    .then(() => {
                        this.proteinFileModel = []; // Reset the file input after upload
                    })
                    .catch(error => {
                        console.error("Upload failed", error);
                    });
            } else {
                console.log("No files selected");
            }
        },
        async handleProteinUpload(files: File[], additionalMetaDataArray?: Record<string, string>[]) {
            this.$q.loading.show({
                spinner: QSpinnerOrbit,
                message: `Uploading protein files`
            });
            const workflowStore = useWorkflowStore();

            try {
                for (let index = 0; index < files.length; index++) {
                    const file = files[index];
                    const metaData = additionalMetaDataArray ? additionalMetaDataArray[index] : undefined;
                    await workflowStore.uploadProteinToExperiment(this.experimentId, '', file, undefined, metaData);
                }
            } catch (error) {
                console.error('Error uploading protein files:', error);
            } finally {
                this.uploadProteinDialog = false;
                this.$q.loading.hide();
            }
        },
        async deleteProtein(proteinToDelete: ProteinResponse) {
            const workflowStore = useWorkflowStore();
            this.$q.loading.show({
                spinner: QSpinnerOrbit,
                message: `Deleting protein`
            });
            try {
                await workflowStore.deleteProteinFromExperiment(proteinToDelete.id);
                this.proteins = workflowStore.proteins;
            } catch (error) {
                console.error('Error deleting protein:', error);
            }
            finally {
                this.$q.loading.hide();
            }
        },
        extractStructureId(url: string): string | null {
            const regex = /\/structure\/(\w+)$/;
            const match = url.match(regex);
            return match ? match[1] : null;
        },
    }
}
)
</script>