<template>
  <q-list bordered separator class="bg-black">
    <q-item-label header>Upload ligands</q-item-label>
    <q-item>
      <q-btn size="md" class="full-width q-pm-sm" push color="info" @click="uploadLigandDialog = true">
        + upload ligands
        <q-tooltip>Add ligands to the experiment</q-tooltip>
      </q-btn>
    </q-item>
    <q-expansion-item
      expand-separator
      v-for="ligand in ligands"
      :key="ligand.ligand_id"
      clickable
      :label="ligand.ligand_name"
      @show="() => getLigandData(ligand)"
    >
      <q-btn color="negative" @click="deleteLigand(ligand)" icon="delete" flat>
        Delete ligand
      </q-btn>
      <LigandDetail v-if="ligand.data" :originalLigand="ligand"> </LigandDetail>
    </q-expansion-item>
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
</template>

<script lang="ts">
import LigandDetail from "src/features/drug_discovery/components/ligands/LigandDetail.vue";
import {useDrugDiscoveryStore} from "src/features/drug_discovery/storage";
import {QSpinnerOrbit} from "quasar";
import {defineComponent} from "vue";
import {ExtendedLigandMetaData, ExtendedTargetMetaData} from "../targets/types";

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
        uploadingLigandFiles: [],
        selectedLigand: null,
      };
    },
  computed: {
    ligands(): ExtendedLigandMetaData[] {
      return this.originalLigands.map(ligand => ({
        ...ligand
      })) as ExtendedLigandMetaData[];
    },
  },
    methods: {
      async handleLigandFileUpload() {
        const store = useDrugDiscoveryStore();

        for (let file of this.uploadingLigandFiles) {
          const ligand = await store.uploadLigandToExperiment(this.experimentId, file);
          const new_ligand = {
            ...ligand,
            jobs: [],
            loadingLigandData: false
          } as ExtendedLigandMetaData
          this.ligands.push(new_ligand);
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
