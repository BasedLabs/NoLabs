<template>
  <q-list bordered separator class="bg-black">
    <q-item-label header>Upload Ligands</q-item-label>
    <q-item>
      <q-btn push  class="full-width" color="info" @click="uploadToAllTargets = true; uploadLigandDialog = true;">
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

      <div v-if="target.loadingLigands">
        <q-spinner color="primary" />
      </div>
      <q-list bordered separator v-else>
        <q-item-label header>Ligands</q-item-label>
        <q-item>
          <q-btn push class="full-width" color="info" @click="openLigandUploadDialog(target)">
            + add ligands
          </q-btn>
        </q-item>
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

          <q-card-section>
            <q-btn
              label="Register a job with current parameters"
              color="info"
              class="q-pa-md"
              @click="() => registerJob(target, ligand)"
            />
          </q-card-section>
          <q-card-section>
            <div class="text-subtitle2 q-pa-md">Docking Jobs:</div>
            <q-list v-if="ligand.jobs && ligand.jobs.length > 0">
              <q-item v-for="(job, jobIndex) in ligand.jobs" :key="job.job_id" class="q-pa-md">
                <q-item-section>Job ID: {{ job.job_id }}</q-item-section>
                <q-item-section side top>
                  <q-btn icon="delete" color="negative" flat @click="deleteJob(target, ligand, job.job_id, jobIndex)" />
                </q-item-section>
              </q-item>
            </q-list>
          </q-card-section>
        </q-expansion-item>
      </q-list>
    </q-expansion-item>
  </q-list>

  <q-dialog v-model="uploadLigandDialog">
    <q-card>
      <q-card-section>
        <div class="text-h6">Upload Ligand Files</div>
        <q-file v-model="uploadingLigandFiles" accept=".sdf" multiple label="Choose files" />
      </q-card-section>

      <q-card-actions align="right">
        <q-btn flat label="Close" color="negative" @click="uploadLigandDialog = false" />
        <q-btn flat label="Upload" color="positive" @click="() => handleLigandFileUpload(selectedTarget!)" />
      </q-card-actions>
    </q-card>
  </q-dialog>
</template>

<script lang="ts">
import TargetDetail from "src/features/drug_discovery/components/targets/TargetDetail.vue";
import LigandDetail from "src/features/drug_discovery/components/ligands/LigandDetail.vue";
import {TargetMetaData} from "src/api/client/models/TargetMetaData.ts";
import {LigandMetaData} from "src/api/client/models/LigandMetaData.ts";
import {ResultMetaData} from "src/api/client/models/ResultMetaData.ts";
import {useDrugDiscoveryStore} from "src/features/drug_discovery/storage";
import {QSpinnerOrbit} from "quasar";
import {defineComponent} from "vue";
import {ExtendedLigandMetaData, ExtendedTargetMetaData} from "./types";

export default defineComponent({
  name: "TargetsLigandsList",
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
      type: Array<TargetMetaData>,
      required: true,
    },
  },
  data() {
    return {
      targets: this.originalTargets.map(target => ({
        ...target,
        loadingLigands: false,
        ligands: [],
        loadingTargetData: false,
      })) as ExtendedTargetMetaData[],
      uploadTargetDialog: false,
      uploadingTargetFiles: [],
      uploadLigandDialog: false,
      uploadingLigandFiles: [] as File[],
      selectedTarget: null as ExtendedTargetMetaData | null,
      uploadToAllTargets: false
    };
  },
  methods: {
    openLigandUploadDialog(target: ExtendedTargetMetaData) {
      this.selectedTarget = target;
      this.uploadLigandDialog = true;
    },
    async selectTarget(target: ExtendedTargetMetaData) {
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: `Selecting target`
      });
      try {
        await this.getLigandsForTarget(target);
        await this.getTargetData(target);
      }
      finally{
        this.$q.loading.hide();
      }
    },
    async getLigandsForTarget(target: ExtendedTargetMetaData) {
      if (target.ligands && target.ligands.length > 0) {
        return;
      }
      const store = useDrugDiscoveryStore();
      target.loadingLigands = true;
      try {
        const originalLigands = await store.fetchLigandsForTarget(
          this.experimentId,
          target.target_id
        );
        target.ligands = originalLigands?.map(ligand => ({
          ...ligand,
          jobs: [],
          loadingLigandData: false,
        })) as ExtendedLigandMetaData[];
      } catch (error) {
        console.error("Error fetching ligands:", error);
      } finally {
        target.loadingLigands = false;
      }
    },
    async handleLigandFileUpload(target: ExtendedTargetMetaData) {
      const store = useDrugDiscoveryStore();
      const targets = this.uploadToAllTargets ? this.targets : [target];

      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: `Uploading ligand`
      });

      try {
        for (let t of targets) {
          for (let file of this.uploadingLigandFiles) {
            const ligand = await store.uploadLigandToTarget(this.experimentId, t.target_id, file);
            const new_ligand = {
              ...ligand,
              jobs: [],
              loadingLigandData: false
            } as ExtendedLigandMetaData
            t.ligands.push(new_ligand);
            await this.registerJob(t, new_ligand);
          }
        }
        this.uploadLigandDialog = false;
        this.uploadToAllTargets = false; // Reset the flag
      }
      finally{
        this.$q.loading.hide();
      }
    },
    async registerJob(target: ExtendedTargetMetaData, ligand: ExtendedLigandMetaData) {
      const store = useDrugDiscoveryStore();
      const response = await store.registerDockingJob(this.experimentId, target.target_id, ligand.ligand_id);
      if (response && response.job_id) {
        ligand.jobs.push({ job_id: response.job_id });
      }
    },
    async deleteJob(target: ExtendedTargetMetaData, ligand: ExtendedLigandMetaData, job_id: string, jobIndex: number) {
      const store = useDrugDiscoveryStore();
      const response = await store.deleteDockingJob(this.experimentId, target.target_id, ligand.ligand_id, job_id);
      if (response) {
        ligand.jobs.splice(jobIndex, 1); // Remove the job from the list
      }
    },
    async getTargetData(target: ExtendedTargetMetaData) {
      const store = useDrugDiscoveryStore();
      target.loadingTargetData = true;
      try {
        const data = await store.fetchTargetData(this.experimentId, target.target_id);
        target.data = {
          proteinSequence: data?.sequence,
          pdbContents: data?.pdbContents,
        };
      } catch (error) {
        console.error('Error fetching target data:', error);
      } finally {
        target.loadingTargetData = false;
      }
    },
    async getLigandData(target: ExtendedTargetMetaData, ligand: ExtendedLigandMetaData) {
      const store = useDrugDiscoveryStore();
      ligand.loadingLigandData = true;
      try {
        const data = await store.fetchLigandData(this.experimentId, target.target_id, ligand.ligand_id);
        ligand.data = {
          sdf_file: data?.ligandSdf,
          smiles: data?.ligandSmiles,
        };
        const jobs = await store.getDockingResultsListForTargetLigand(this.experimentId, target.target_id, ligand.ligand_id);
        ligand.jobs = jobs?.results_list as ResultMetaData[];// Store the jobs in the ligand object
      } catch (error) {
        console.error('Error fetching ligand data:', error);
      } finally {
        ligand.loadingLigandData = false;
      }
    },
    async deleteTarget(targetToDelete: TargetMetaData) {
      const store = useDrugDiscoveryStore();
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: `Loading ligands`
      });
      try {
        await store.deleteTargetFromExperiment(this.experimentId, targetToDelete.target_id);
        this.targets = this.targets.filter(target => target.target_id !== targetToDelete.target_id);
      } catch (error) {
        console.error('Error deleting target:', error);
      }
      finally {
        this.$q.loading.hide();
      }
    },
    async deleteLigand(target: TargetMetaData, ligandToDelete: LigandMetaData) {
      const store = useDrugDiscoveryStore();
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: `Deleting ligand`
      });
      try {
        await store.deleteLigandFromTarget(this.experimentId, target.target_id, ligandToDelete.ligand_id);
        target.ligands = target.ligands.filter(ligand => ligand.ligand_id !== ligandToDelete.ligand_id);
      } catch (error) {
        console.error('Error deleting ligand:', error);
      }
      finally{
        this.$q.loading.hide();
      }
    },
  },
}
)
</script>
