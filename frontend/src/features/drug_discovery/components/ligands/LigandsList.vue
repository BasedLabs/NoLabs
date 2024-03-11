<template>

  <q-list bordered separator style="width: 100%;" class="bg-black rounded-borders">
    <q-item-label header class="text-white">Upload ligands</q-item-label>
    <q-item>
      <q-btn size="md" class="full-width q-pm-sm" push color="info" @click="uploadLigandDialog = true">
        + upload ligands
        <q-tooltip>Add ligands to the experiment</q-tooltip>
      </q-btn>
    </q-item>
    <q-item v-for="ligand in ligands" :key="ligand.ligand_id" clickable v-ripple class="q-mb-sm" @click="showLigandDetailDialog(ligand)">
      <q-card class="q-pa-md full-width">
        <q-card-section>
          <div>{{ ligand.ligand_name }}</div>
          <div class="q-pt-sm q-pb-sm"><a class="text-light-blue" :href="ligand.link" target="_blank">{{ ligand.link }}</a></div>
          <q-item-label v-if="ligand.data" caption>SMILES: {{ ligand.data.smiles }}</q-item-label>
        </q-card-section>
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
        <q-file v-model="uploadingLigandFiles" accept=".sdf" multiple label="Choose files" />
      </q-card-section>

      <q-card-actions align="right">
        <q-btn flat label="Close" color="negative" @click="uploadLigandDialog = false" />
        <q-btn flat label="Upload" color="positive" @click="() => handleLigandFileUpload()" />
      </q-card-actions>
    </q-card>
  </q-dialog>

  <q-dialog v-model="ligandDetailDialogVisible">
    <q-card>
      <q-card-section>
        <LigandDetail :originalLigand="selectedLigand" v-if="selectedLigand && selectedLigand.data"></LigandDetail>
      </q-card-section>
      <q-card-actions align="right">
        <q-btn flat label="Close" color="negative" @click="ligandDetailDialogVisible = false" />
      </q-card-actions>
    </q-card>
  </q-dialog>

</template>

<script lang="ts">
import LigandDetail from "src/features/drug_discovery/components/ligands/LigandDetail.vue";
import {useDrugDiscoveryStore} from "src/features/drug_discovery/storage";
import {QSpinnerOrbit} from "quasar";
import {defineComponent} from "vue";
import {ExtendedLigandMetaData} from "../targets/types";
import {useBioBuddyStore} from "../../../biobuddy/storage";
import {ChemBLData, FunctionCallReturnData} from "../../../../api/client";

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
      originalLigands: {
        type: Array as () => ExtendedLigandMetaData[],
        required: true,
      },
    },
    data() {
      return {
        uploadLigandDialog: false,
        uploadingLigandFiles: [] as File[],
        selectedLigand: null as ExtendedLigandMetaData | null,
        chemblQueryCallback:  null as ((data: FunctionCallReturnData) => void) | null,
        ligandDetailDialogVisible: false,
      };
    },
    computed: {
      ligands(): ExtendedLigandMetaData[] {
        return this.originalLigands.map(ligand => ({
          ...ligand
        })) as ExtendedLigandMetaData[];
      },
    },
    mounted() {
      const bioBuddyStore = useBioBuddyStore();
      this.chemblQueryCallback = (data: FunctionCallReturnData) => {
        const metaDatasArray: Record<string, string>[] = [];
        for (let dataObject of data.files){
          const chemBLObject = dataObject as ChemBLData;
          const file = new File([new Blob([chemBLObject.content!])], chemBLObject.metadata.chembl_id + ".sdf")

          this.uploadingLigandFiles.push(file);
          metaDatasArray.push(chemBLObject.metadata);
        }
        this.handleLigandFileUpload(metaDatasArray);
      };
      bioBuddyStore.addQueryChemblEventHandler(this.chemblQueryCallback);

      for (let ligand of this.ligands) {
        this.getLigandData(ligand);
      }
    },
    methods: {
      showLigandDetailDialog(ligand: ExtendedLigandMetaData) {
        this.selectedLigand = ligand;
        this.ligandDetailDialogVisible = true;
      },
      async handleLigandFileUpload(additionalMetaDataArray?: Record<string, string>[]) {
        const store = useDrugDiscoveryStore();

        for (let index = 0; index < this.uploadingLigandFiles.length; index++) {
          const file = this.uploadingLigandFiles[index];
          const metaData = additionalMetaDataArray ? additionalMetaDataArray[index] : undefined;
          await store.uploadLigandToExperiment(this.experimentId, file, metaData);
        }
        this.uploadLigandDialog = false;
      },
      async getLigandData(ligand: ExtendedLigandMetaData) {
        const store = useDrugDiscoveryStore();
        ligand.loadingLigandData = true;
        try {
          await store.fetchLigandDataForExperiment(this.experimentId, ligand.ligand_id);
        } catch (error) {
          console.error('Error fetching ligand data:', error);
        } finally {
          ligand.loadingLigandData = false;
        }
      },
      async deleteLigand(ligandToDelete: ExtendedLigandMetaData) {
        const store = useDrugDiscoveryStore();
        this.$q.loading.show({
          spinner: QSpinnerOrbit,
          message: `Deleting ligand`
        });
        try {
          await store.deleteLigandFromExperiment(this.experimentId, ligandToDelete.ligand_id);
          this.ligands = this.ligands.filter(ligand => ligand.ligand_id !== ligandToDelete.ligand_id);
        } catch (error) {
          console.error('Error deleting ligand:', error);
        } finally {
          this.$q.loading.hide();
        }
      },
    }
  }
)
</script>
