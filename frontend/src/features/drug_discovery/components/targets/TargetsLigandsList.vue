<template>
  <q-list bordered separator>
    <q-item-label header>Upload Ligands</q-item-label>
    <q-item>
      <q-btn push  class="full-width" color="green" @click="uploadToAllTargets = true; uploadLigandDialog = true;">
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
      <q-btn color="negative" @click="deleteTarget(target)" icon="delete" flat>
        Delete target
      </q-btn>
      <q-expansion-item
        popup
        expand-separator
        :label="'Show target data'"
        :content-inset-level="1"
      ><TargetDetail v-if="target.data" :experiment-id="experimentId" :original-target="target"> </TargetDetail>
      </q-expansion-item>

      <q-btn push class="full-width"  color="green" @click="openLigandUploadDialog(target)">
        + add ligands
      </q-btn>
      <div v-if="target.loadingLigands">
        <q-spinner color="primary" />
      </div>
      <q-list bordered separator v-else>
        <q-item-label header>Ligands</q-item-label>
        <q-expansion-item
          expand-separator
          v-for="ligand in target.ligands"
          :key="ligand.ligand_id"
          clickable
          :label="ligand.ligand_name"
          @show="() => getLigandData(target, ligand)"
        >
          <q-btn color="negative" @click="deleteLigand(target, ligand)" icon="delete" flat>
            Delete ligand
          </q-btn>
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
        <q-file v-model="uploadingTargetFiles" accept=".fasta" multiple label="Choose files" />
      </q-card-section>

      <q-card-actions align="right">
        <q-btn flat label="Close" color="negative" @click="uploadTargetDialog = false" />
        <q-btn flat label="Upload" color="positive" @click="handleTargetFileUpload" />
      </q-card-actions>
    </q-card>
  </q-dialog>

  <q-dialog v-model="uploadLigandDialog">
    <q-card>
      <q-card-section>
        <div class="text-h6">Upload Ligand Files</div>
        <q-file v-model="uploadingLigandFiles" accept=".sdf" multiple label="Choose files" />
      </q-card-section>

      <q-card-actions align="right">
        <q-btn flat label="Close" color="primary" @click="uploadLigandDialog = false" />
        <q-btn flat label="Upload" color="primary" @click="() => handleLigandFileUpload(selectedTarget)" />
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
      uploadLigandDialog: false,
      uploadingLigandFiles: [],
      selectedTarget: null,
      uploadToAllTargets: false,
    };
  },
  methods: {
    openLigandUploadDialog(target) {
      this.selectedTarget = target;
      this.uploadLigandDialog = true;
    },
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
    async handleLigandFileUpload(target) {
      const store = useDrugDiscoveryStore();
      const targets = this.uploadToAllTargets ? this.targets : [target];

      for (let t of targets) {
        for (let file of this.uploadingLigandFiles) {
          const ligand = await store.uploadLigandToTarget(this.experimentId, t.target_id, file);
          t.ligands.push(ligand);
        }
      }
      this.uploadLigandDialog = false;
      this.uploadToAllTargets = false; // Reset the flag
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
    async deleteLigand(target, ligandToDelete) {
      const store = useDrugDiscoveryStore();
      try {
        await store.deleteLigandFromTarget(this.experimentId, target.target_id, ligandToDelete.ligand_id);
        target.ligands = target.ligands.filter(ligand => ligand.ligand_id !== ligandToDelete.ligand_id);
      } catch (error) {
        console.error('Error deleting ligand:', error);
      }
    },
  },
};
</script>
