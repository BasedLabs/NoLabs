<template>
  <q-list bordered style="{ width: scrollContentWidth + 'px' }" class="bg-black rounded-borders">
    <q-item-label header class="text-white text-h6">Ligands List</q-item-label>
    <q-item>
      <q-btn size="md" class="full-width q-pm-sm" push color="info" @click="uploadLigandDialog = true">
        + upload ligands
        <q-tooltip>Add ligands to the experiment</q-tooltip>
      </q-btn>
    </q-item>
    <div :style="{ width: scrollContentWidth + 'px' }">
      <q-scroll-area v-if="filteredLigands.length > 0" visible :thumbStyle="thumbStyle" :barStyle="barStyle" style="height: 60vh;">
        <div ref="scrollContent" class="scroll-content">
          <q-item v-for="ligand in filteredLigands" :key="ligand.id" clickable v-ripple class="q-mb-sm"
            @click="showLigandDetailDialog(ligand)">
            <q-card class="q-pa-md full-width flex-row justify-between" bordered>
              <div class="row">
                <div class="col-4">
                  <img v-if="ligand.image" class="rounded-borders" :width="100" :height="100"
                    :src="'data:image/png;base64,' + ligand.image" alt="Ligand Structure" />
                </div>
                <div class="col-8">
                  <q-card-section>
                    <div>{{ ligand.name }}</div>
                    <q-item-label v-if="ligand.smiles_content" caption>SMILES: {{
                      ligand.smiles_content
                    }}
                    </q-item-label>
                  </q-card-section>
                  <q-card-section class="flex-row">
                    <q-btn icon="delete" flat @click.stop="deleteLigand(ligand)" color="negative">
                      <q-tooltip>Delete ligand</q-tooltip>
                    </q-btn>
                    <q-btn icon="file_download" flat @click.stop="downloadLigand(ligand)" color="info">
                      <q-tooltip>Download</q-tooltip>
                    </q-btn>
                  </q-card-section>
                </div>
              </div>
            </q-card>
          </q-item>
        </div>
      </q-scroll-area>
    </div>
  </q-list>

  <q-dialog v-model="uploadLigandDialog">
    <q-card>
      <q-card-section>
        <div class="text-h6">Upload Ligand Files</div>
        <q-file v-model="fileModel" accept=".sdf" multiple label="Choose files" />
      </q-card-section>

      <q-card-actions align="right">
        <q-btn flat label="Close" color="negative" @click="uploadLigandDialog = false" />
        <q-btn flat label="Upload" color="positive" @click="uploadFiles" />
      </q-card-actions>
    </q-card>
  </q-dialog>

  <q-dialog v-model="ligandDetailDialogVisible">
    <q-card>
      <q-card-section>
        <LigandDetail :ligandId="selectedLigand.id" v-if="selectedLigand && selectedLigand.sdf_content">
        </LigandDetail>
      </q-card-section>
      <q-card-actions align="right">
        <q-btn flat label="Close" color="negative" @click="ligandDetailDialogVisible = false" />
      </q-card-actions>
    </q-card>
  </q-dialog>
</template>


<script lang="ts">
import LigandDetail from "src/features/workflow/components/nodeTemplates/dataSourceNodes/ligands/LigandDetail.vue";
import { useWorkflowStore } from "src/features/workflow/components/storage";
import { QSpinnerOrbit } from "quasar";
import { defineComponent, onMounted, ref } from "vue";
import { LigandContentResponse } from "src/refinedApi/client";
import { useBioBuddyStore } from "src/features/biobuddy/storage";
import { ChemBLData } from "src/refinedApi/client";

export default defineComponent({
  name: "LigandsList",
  components: {
    LigandDetail,
  },
  props: {
    experimentId: {
      type: String,
      required: true,
    },
    nodeId: {
      type: String,
      required: true,
    }
  },
  data() {
    return {
      uploadLigandDialog: false,
      selectedLigand: null as LigandContentResponse | null,
      fileModel: [],
      chemblQueryCallback: null as ((data: { files: ChemBLData[] }) => void) | null,
      ligandDetailDialogVisible: false,
      ligands: [] as LigandContentResponse[],
      thumbStyle: {
        right: '4px',
        borderRadius: '7px',
        backgroundColor: '#027be3',
        width: '4px',
        opacity: 0.75
      },
      barStyle: {
        right: '2px',
        borderRadius: '9px',
        backgroundColor: '#027be3',
        width: '8px',
        opacity: 0.2
      },
      scrollContentWidth: 0
    };
  },
  async mounted() {
    const workflowStore = useWorkflowStore();
    await workflowStore.getAllLigands(this.experimentId);
    this.ligands = workflowStore.ligands;

    const bioBuddyStore = useBioBuddyStore();
    this.chemblQueryCallback = async (data: { files: ChemBLData[] }) => {
      const files: File[] = [];
      const metaDatasArray: Record<string, string>[] = [];

      for (let chemBLObject of data.files) {
        const file = new File([new Blob([chemBLObject.content!])], chemBLObject.metadata.pref_name + ".sdf");
        files.push(file);
        metaDatasArray.push(chemBLObject.metadata);
      }
      await this.handleLigandFileUpload(files, metaDatasArray);
    };
    bioBuddyStore.addQueryChemblEventHandler(this.chemblQueryCallback);
    this.updateScrollContentWidth();
  },
  computed: {
    filteredLigands() {
      const workflowStore = useWorkflowStore();
      const nodeData = workflowStore.getNodeById(this.nodeId);
      const ligandIds = nodeData?.data.defaults[0]?.value || [];
      return this.ligands.filter(ligand => ligandIds.includes(ligand.id));
    }
  },
  methods: {
    showLigandDetailDialog(ligand: LigandContentResponse) {
      this.selectedLigand = ligand;
      this.ligandDetailDialogVisible = true;
      this.updateScrollContentWidth();
    },
    uploadFiles() {
      if (this.fileModel.length > 0) {
        this.handleLigandFileUpload(this.fileModel)
          .then(() => {
            this.fileModel = []; // Reset/clear the file input after upload
          })
          .catch(error => {
            console.error("Upload failed", error);
          });
      } else {
        console.log("No files selected");
      }
    },
    async handleLigandFileUpload(files: File[], additionalMetaDataArray?: Record<string, string>[]) {
      const workflowStore = useWorkflowStore();

      for (let index = 0; index < files.length; index++) {
        const file = files[index];
        const metaData = additionalMetaDataArray ? additionalMetaDataArray[index] : undefined;
        const newLigand = await workflowStore.uploadLigandToExperiment(this.experimentId, this.nodeId, file.name, undefined, file, metaData);
        this.ligands = workflowStore.ligands;
      }
      this.uploadLigandDialog = false;
      this.updateScrollContentWidth();
    },
    async deleteLigand(ligandToDelete: LigandContentResponse) {
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: `Deleting ligand`
      });
      const workflowStore = useWorkflowStore();
      try {
        await workflowStore.deleteLigandFromExperiment(this.nodeId, ligandToDelete.id);
        this.ligands = workflowStore.ligands;
      } catch (error) {
        console.error('Error deleting ligand:', error);
      } finally {
        this.$q.loading.hide();
      }
    },
    downloadLigand(ligand: LigandContentResponse) {
      if (!ligand.sdf_content) {
        console.error('SDF file content is missing for this ligand.');
        return;
      }

      const sdfContent = ligand.sdf_content;
      const blob = new Blob([sdfContent], { type: 'chemical/x-mdl-sdfile' });
      const url = window.URL.createObjectURL(blob);

      const a = document.createElement('a');
      a.href = url;
      a.download = `${ligand.name}.sdf`;
      document.body.appendChild(a);
      a.click();
      window.URL.revokeObjectURL(url);
      a.remove();
    },
    updateScrollContentWidth() {
      this.$nextTick(() => {
        const scrollContent = this.$refs.scrollContent as HTMLElement;
        this.scrollContentWidth = scrollContent ? scrollContent.scrollWidth : 0;
      });
    }
  },
  unmounted() {
    const bioBuddyStore = useBioBuddyStore();
    const index = bioBuddyStore.queryChemblEventHandlers.indexOf(this.chemblQueryCallback!);
    if (index !== -1) {
      bioBuddyStore.queryChemblEventHandlers.splice(index, 1);
    }
  }
});
</script>


<style scoped>
.scroll-content {
  display: flex;
  flex-direction: column;
  width: max-content;
}

.q-list {
  overflow-x: auto;
}
</style>
