<template>
  <q-list bordered separator>
    <q-item-label header>Upload targets</q-item-label>
    <q-expansion-item
      popup
      expand-separator
      :content-inset-level="1"
      v-for="target in targets"
      :key="target.target_id"
      :label="target.target_name"
      @show="() => selectTarget(target)"
    >
      <q-btn color="negative" @click="deleteTarget(target)" icon="delete" flat>
      Delete target
      </q-btn>
      <TargetDetail v-if="target.data" :experiment-id="experimentId" :original-target="target"> </TargetDetail>
    </q-expansion-item>
  </q-list>
  <q-page-sticky position="bottom-right" :offset="[100, 50]">
    <q-btn size="lg" round color="green" icon="add" @click="uploadTargetDialog = true">
      <q-tooltip>Add targets to the experiment</q-tooltip>
    </q-btn>
  </q-page-sticky>

  <q-dialog v-model="uploadTargetDialog">
    <q-card>
      <q-card-section>
        <div class="text-h6">Upload Target Files</div>
        <q-file v-model="uploadingTargetFiles" accept=".fasta" multiple label="Choose files" />
      </q-card-section>

      <q-card-actions align="right">
        <q-btn flat label="Close" color="negative" @click="uploadTargetDialog = false" />
        <q-btn flat label="Upload" color="positive" @click="handleTargetFileUpload" />
      </q-card-actions>
    </q-card>
  </q-dialog>
</template>

<script>
import TargetDetail from "src/features/drug_discovery/components/targets/TargetDetail.vue";
import { useDrugDiscoveryStore } from "src/features/drug_discovery/storage";

export default {
  name: "TargetsList",
  components: {
    TargetDetail,
  },
  props: {
    experimentId: {
      type: String,
      required: true,
    },
    originalTargets: {
      type: Array,
      required: true,
    },
  },
  data() {
    return {
      targets: this.originalTargets,
      uploadTargetDialog: false,
      uploadingTargetFiles: [],
      uploadLigandDialog: false,
      uploadingLigandFiles: [],
      selectedTarget: null,
      uploadToAllTargets: false,
    };
  },
  methods: {
    async selectTarget(target) {
      await this.getLigandsForTarget(target);
      await this.getTargetData(target);
    },
    async getLigandsForTarget(target) {
      if (target.ligands && target.ligands.length > 0) {
        return;
      }
      const store = useDrugDiscoveryStore();
      target.loadingLigands = true;
      try {
        target.ligands = await store.fetchLigandsForTarget(
          this.experimentId,
          target.target_id
        );
      } catch (error) {
        console.error("Error fetching ligands:", error);
      } finally {
        target.loadingLigands = false;
      }
    },
    async handleTargetFileUpload() {
      const store = useDrugDiscoveryStore();
      for (let file of this.uploadingTargetFiles) {
        await store.uploadTargetToExperiment(this.experimentId, file);
      }
      this.uploadTargetDialog = false;
    },
    async getTargetData(target) {
      const store = useDrugDiscoveryStore();
      target.loadingTargetData = true;
      try {
        const data = await store.fetchTargetData(this.experimentId, target.target_id);
        target.data = {
          proteinSequence: data.sequence,
          pdbContents: data.pdbContents,
        };
      } catch (error) {
        console.error('Error fetching target data:', error);
      } finally {
        target.loadingTargetData = false;
      }
    },
    async deleteTarget(targetToDelete) {
      const store = useDrugDiscoveryStore();
      try {
        await store.deleteTargetFromExperiment(this.experimentId, targetToDelete.target_id);
        this.targets = this.targets.filter(target => target.target_id !== targetToDelete.target_id);
      } catch (error) {
        console.error('Error deleting target:', error);
      }
    },
  },
};
</script>
