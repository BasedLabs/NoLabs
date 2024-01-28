<template>
  <q-list bordered separator>
    <q-item-label header>Targets</q-item-label>
    <q-item>
      <q-btn push class="full-width" color="green">
        + add ligands to all targets
      </q-btn>
    </q-item>
    <q-expansion-item
      popup
      expand-separator
      :content-inset-level="1"
      v-for="target in targets"
      :key="target.target_id"
      :label="target.target_name"
      @show="() => selectTarget(target)"
    >
      <TargetDetail v-if="target.data" :experiment-id="experimentId" :original-target="target"> </TargetDetail>
      <q-btn push class="full-width" color="green"> + add ligand </q-btn>
      <div v-if="target.loadingLigands">
        <q-spinner color="primary" />
      </div>
      <q-list bordered v-else>
        <q-expansion-item
          v-for="ligand in target.ligands"
          :key="ligand.ligand_id"
          clickable
          :label="ligand.ligand_name"
          @show="() => getLigandData(target, ligand)"
        >
          <LigandDetail v-if="ligand.data" :originalLigand="ligand"> </LigandDetail>
        </q-expansion-item>
      </q-list>
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
        <q-file v-model="uploadingTargetFiles" multiple label="Choose files" />
      </q-card-section>

      <q-card-actions align="right">
        <q-btn flat label="Close" color="primary" @click="uploadTargetDialog = false" />
        <q-btn flat label="Upload" color="primary" @click="handleTargetFileUpload" />
      </q-card-actions>
    </q-card>
  </q-dialog>
</template>

<script>
import TargetDetail from "src/features/drug_discovery/components/targets/TargetDetail.vue";
import LigandDetail from "src/features/drug_discovery/components/ligands/LigandDetail.vue";
import { useDrugDiscoveryStore } from "src/features/drug_discovery/storage";

export default {
  name: "TargetsList",
  components: {
    TargetDetail,
    LigandDetail,
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
      this.uploadTargetDialog = false; // Close the dialog after upload
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
    async getLigandData(target, ligand) {
      const store = useDrugDiscoveryStore();
      ligand.loadingLigandData = true;
      try {
        const data = await store.fetchLigandData(this.experimentId, target.target_id, ligand.ligand_id);
        ligand.data = {
          sdf_file: data.ligandSdf,
          smiles: data.ligandSmiles,
        };
      } catch (error) {
        console.error('Error fetching target data:', error);
      } finally {
        ligand.loadingLigandData = false;
      }
    }
  },
};
</script>
