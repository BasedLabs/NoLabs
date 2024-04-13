<template>

  <q-list bordered style="width: 100%;" class="bg-black rounded-borders">
    <q-item-label header class="text-white">Upload ligands</q-item-label>
    <q-item>
      <q-btn size="md" class="full-width q-pm-sm" push color="info" @click="uploadLigandDialog = true">
        + upload ligands
        <q-tooltip>Add ligands to the experiment</q-tooltip>
      </q-btn>
    </q-item>
    <q-item v-for="ligand in ligands" :key="ligand.ligand_id" clickable v-ripple class="q-mb-sm"
            @click="showLigandDetailDialog(ligand)">
      <q-card class="q-pa-md full-width flex-row justify-between" bordered>
        <div class="row">
          <div class="col-4">
            <img v-if="ligand.image" class="rounded-borders" :width="200" :height="200"
                 :src="'data:image/png;base64,' + ligand.image" alt="Ligand Structure"/>
          </div>
          <div class="col-8">
            <q-card-section>
              <div>{{ ligand.ligand_name }}</div>
              <div class="q-pt-sm q-pb-sm"><a class="text-light-blue" :href="ligand.link" target="_blank">{{
                  ligand.link
                }}</a></div>
              <q-item-label v-if="ligand.data" caption>SMILES: {{ ligand.data.smiles }}</q-item-label>
            </q-card-section>
            <q-card-section class="flex-row">
              <q-btn icon="delete" flat @click.stop="deleteLigand(ligand)" color="negative">
                <q-tooltip>Delete ligand</q-tooltip>
              </q-btn>
              <q-btn icon="add" flat @click.stop="openTargetSelectionDialog(ligand)" color="positive">
                <q-tooltip>Attach to targets</q-tooltip>
              </q-btn>
              <q-btn icon="file_download" flat @click.stop="downloadLigand(ligand)" color="info">
                <q-tooltip>Download</q-tooltip>
              </q-btn>
            </q-card-section>
          </div>
        </div>
      </q-card>
    </q-item>
  </q-list>


  <q-page-sticky position="bottom-right" :offset="[100, 50]">
    <q-btn size="lg" push round color="info" icon="add" @click="uploadLigandDialog = true">
      <q-tooltip>Add ligands to the experiment</q-tooltip>
    </q-btn>
  </q-page-sticky>
  <q-dialog v-model="uploadLigandDialog">
    <q-card>
      <q-card-section>
        <div class="text-h6">Upload Ligand Files</div>
        <q-file v-model="fileModel" accept=".sdf" multiple label="Choose files"/>
      </q-card-section>

      <q-card-actions align="right">
        <q-btn flat label="Close" color="negative" @click="uploadLigandDialog = false"/>
        <q-btn flat label="Upload" color="positive" @click="uploadFiles"/>
      </q-card-actions>
    </q-card>
  </q-dialog>

  <q-dialog v-model="ligandDetailDialogVisible">
    <q-card>
      <q-card-section>
        <LigandDetail :originalLigand="selectedLigand" v-if="selectedLigand && selectedLigand.data"></LigandDetail>
      </q-card-section>
      <q-card-actions align="right">
        <q-btn flat label="Close" color="negative" @click="ligandDetailDialogVisible = false"/>
      </q-card-actions>
    </q-card>
  </q-dialog>

  <q-dialog v-model="targetSelectionDialogVisible">
    <q-card>
      <q-card-section>
        <div class="text-h6">Select Targets</div>
        <q-checkbox label="Select All" v-model="allTargetsSelected"></q-checkbox>
        <q-list>
          <q-item v-for="target in targets" :key="target.target_id">
            <q-checkbox v-model="this.selectedTargets" :label="target.target_name" :val="target.target_id"></q-checkbox>
          </q-item>
        </q-list>
      </q-card-section>
      <q-card-actions align="right">
        <q-btn flat label="Add ligands to selected targets" color="positive" @click="addLigandsToSelectedTargets"/>
        <q-btn flat label="Close" color="negative" @click="targetSelectionDialogVisible = false"/>
      </q-card-actions>
    </q-card>
  </q-dialog>

</template>

<script lang="ts">
import LigandDetail from "src/features/drug_discovery/components/ligands/LigandDetail.vue";
import {useDrugDiscoveryStore} from "src/features/drug_discovery/storage";
import {QSpinnerOrbit} from "quasar";
import {defineComponent} from "vue";
import {ExtendedLigandMetaData, ExtendedTargetMetaData} from "../targets/types";
import {useBioBuddyStore} from "../../../biobuddy/storage";
import {ChemBLData} from "src/api/client";

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
      },
      data() {
        return {
          uploadLigandDialog: false,
          selectedLigand: null as ExtendedLigandMetaData | null,
          fileModel: [],
          chemblQueryCallback: null as ((data: { files: ChemBLData[] }) => void) | null,
          ligandDetailDialogVisible: false,
          targetSelectionDialogVisible: false,
          selectedTargets: [] as string[],
          selectAllTargets: false,
          ligands: [] as ExtendedLigandMetaData[],
          targets: [] as ExtendedTargetMetaData[],
        };
      },
      computed: {
        drugDiscoveryStore() {
          return useDrugDiscoveryStore();
        },
        allTargetsSelected: {
          get() {
            return this.selectedTargets.length === this.drugDiscoveryStore.targets.length &&
                this.selectedTargets.length > 0;
          },
          set(value: boolean) {
            this.toggleSelectAll(value);
          }
        }
      },
      async mounted() {
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

        this.ligands = await this.drugDiscoveryStore.fetchLigandsForExperiment(this.experimentId) as ExtendedLigandMetaData[];

        for (let ligand of this.ligands) {
          await this.getLigandData(ligand);
        }
      },
      methods: {
        showLigandDetailDialog(ligand: ExtendedLigandMetaData) {
          this.selectedLigand = ligand;
          this.ligandDetailDialogVisible = true;
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
          for (let index = 0; index < files.length; index++) {
            const file = files[index];
            const metaData = additionalMetaDataArray ? additionalMetaDataArray[index] : undefined;
            const newLigand = await this.drugDiscoveryStore.uploadLigandToExperiment(this.experimentId, file, metaData);
            this.ligands.push(newLigand!);
          }
          this.uploadLigandDialog = false;
        },
        async getLigandData(ligand: ExtendedLigandMetaData) {
          ligand.loadingLigandData = true;
          try {
            ligand.data = await this.drugDiscoveryStore.fetchLigandDataForExperiment(this.experimentId, ligand.ligand_id);
          } catch (error) {
            console.error('Error fetching ligand data:', error);
          } finally {
            ligand.loadingLigandData = false;
          }
        },
        async deleteLigand(ligandToDelete: ExtendedLigandMetaData) {
          this.$q.loading.show({
            spinner: QSpinnerOrbit,
            message: `Deleting ligand`
          });
          try {
            await this.drugDiscoveryStore.deleteLigandFromExperiment(this.experimentId, ligandToDelete.ligand_id);
            this.ligands = this.ligands.filter(ligand => ligand.ligand_id !== ligandToDelete.ligand_id);
          } catch (error) {
            console.error('Error deleting ligand:', error);
          } finally {
            this.$q.loading.hide();
          }
        },
        toggleSelectAll(value: boolean) {
          if (value) {
            // Replace the array with a new one containing all IDs to ensure reactivity
            this.selectedTargets = this.drugDiscoveryStore.targets.map(target => target.target_id);
          } else {
            // Clear the array by replacing it with a new empty array
            this.selectedTargets = [];
          }
        },

        async openTargetSelectionDialog(ligand: ExtendedLigandMetaData) {
          this.selectedLigand = ligand;
          this.targetSelectionDialogVisible = true;
          this.selectedTargets = [];
          this.selectAllTargets = false;
          this.targets = await this.drugDiscoveryStore.fetchTargetsForExperiment(this.experimentId) as ExtendedTargetMetaData[];
        },

        async addLigandsToSelectedTargets() {
          const store = useDrugDiscoveryStore();

          if (this.selectedLigand && this.selectedLigand.data && this.selectedLigand.data.sdf_file) {
            for (let target_id of this.selectedTargets) {
              const sdfFileContent = this.selectedLigand.data.sdf_file;

              const file = new File([new Blob([sdfFileContent])], this.selectedLigand.ligand_name + ".sdf")

              try {
                await store.uploadLigandToTarget(this.experimentId, target_id, file);
              } catch (error) {
                console.error('Error uploading ligand to target:', error);
              }
            }
          }

          this.targetSelectionDialogVisible = false;
          this.selectedTargets = [];
          this.selectAllTargets = false;
        },
        downloadLigand(ligand: ExtendedLigandMetaData) {
          if (!ligand.data || !ligand.data.sdf_file) {
            console.error('SDF file content is missing for this ligand.');
            return;
          }

          const sdfContent = ligand.data.sdf_file;
          const blob = new Blob([sdfContent], {type: 'chemical/x-mdl-sdfile'});
          const url = window.URL.createObjectURL(blob);

          const a = document.createElement('a');
          a.href = url;
          a.download = `${ligand.ligand_name}.sdf`;
          document.body.appendChild(a);
          a.click();
          window.URL.revokeObjectURL(url);
          a.remove();
        },
      },
      unmounted() {
        const bioBuddyStore = useBioBuddyStore();
        const index = bioBuddyStore.queryChemblEventHandlers.indexOf(this.chemblQueryCallback!);
        if (index !== -1) {
          bioBuddyStore.queryChemblEventHandlers.splice(index, 1);
        }
      }
    }
)
</script>
