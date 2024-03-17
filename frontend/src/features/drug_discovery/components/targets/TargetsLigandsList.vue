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
          @show="() => getLigandDataForTarget(target, ligand)"
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
        <!-- Lone Ligands Selection -->
        <div class="q-mt-md">
          <q-card-actions align="right">
            <q-btn flat label="Close" color="negative" @click="uploadLigandDialog = false" />
            <q-btn flat label="Upload" color="positive" @click="() => handleLigandFileUpload(selectedTarget!)" />
          </q-card-actions>
          <div>Or select some of the uploaded ligands:</div>
          <q-list>
            <q-item v-for="ligand in loneLigands" :key="ligand.ligand_id">
              <q-checkbox v-model="selectedLoneLigands" :label="ligand.ligand_name" :val="ligand"></q-checkbox>
            </q-item>
          </q-list>
          <q-card-actions align="right">
            <q-btn flat label="Close" color="negative" @click="uploadLigandDialog = false" />
            <q-btn flat label="Upload" color="positive" @click="() => handleLigandFileUpload(selectedTarget!)" />
          </q-card-actions>
        </div>
      </q-card-section>
    </q-card>
  </q-dialog>

</template>

<script lang="ts">
import TargetDetail from "src/features/drug_discovery/components/targets/TargetDetail.vue";
import LigandDetail from "src/features/drug_discovery/components/ligands/LigandDetail.vue";
import {TargetMetaData} from "src/api/client/models/TargetMetaData.ts";
import {useDrugDiscoveryStore} from "src/features/drug_discovery/storage";
import {QSpinnerOrbit} from "quasar";
import {defineComponent} from "vue";
import {ExtendedLigandMetaData, ExtendedTargetMetaData} from "./types";
import {JobMetaData} from "../../../../api/client";

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
      uploadTargetDialog: false,
      uploadingTargetFiles: [],
      uploadLigandDialog: false,
      uploadingLigandFiles: [] as File[],
      selectedTarget: null as ExtendedTargetMetaData | null,
      uploadToAllTargets: false,
      selectedLoneLigands: [] as ExtendedLigandMetaData[],
      targets: [] as ExtendedTargetMetaData[],
      loneLigands: [] as ExtendedLigandMetaData[]
    };
  },
  async mounted() {
    this.targets = await this.drugDiscoveryStore.fetchTargetsForExperiment(this.experimentId) as ExtendedTargetMetaData[];
    this.loneLigands = await this.drugDiscoveryStore.fetchLigandsForExperiment(this.experimentId) as ExtendedLigandMetaData[];
  },
  computed: {
      drugDiscoveryStore() {
        return useDrugDiscoveryStore();
      },
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
        target.ligands = await this.getLigandsForTarget(target) as ExtendedLigandMetaData[];
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
      target.loadingLigands = true;
      try {
        const originalLigands = await this.drugDiscoveryStore.fetchLigandsForTarget(
          this.experimentId,
          target.target_id
        );
        return originalLigands?.map(ligand => ({
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
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: `Uploading ligands...`
      });

      try {
        // Handle file uploads
        const targets = this.uploadToAllTargets ? this.targets : [target];
        for (let t of targets) {
          for (let file of this.uploadingLigandFiles) {
            const ligandMetaData = await this.drugDiscoveryStore.uploadLigandToTarget(this.experimentId, t.target_id, file);
            t.ligands.push(ligandMetaData!);
          }
        }

        // Handle selected lone ligands
        for (let ligand of this.selectedLoneLigands) {
          if (!ligand.data?.sdf_file) {
            ligand.data = await this.drugDiscoveryStore.fetchLigandDataForExperiment(this.experimentId, ligand.ligand_id);
          }
          const file = new File([new Blob([ligand.data?.sdf_file!])], `${ligand.ligand_name}.sdf`);
          for (let t of targets) {
            const ligandMetaData = await this.drugDiscoveryStore.uploadLigandToTarget(this.experimentId, t.target_id, file);
            t.ligands.push(ligandMetaData!);
          }
        }
      } catch (error) {
        console.error('Error during ligand upload:', error);
      } finally {
        this.$q.loading.hide();
        this.uploadLigandDialog = false;
        this.uploadToAllTargets = false; // Reset the flag
        this.selectedLoneLigands = []; // Clear selected lone ligands
      }
    },
    async registerJob(target: ExtendedTargetMetaData, ligand: ExtendedLigandMetaData) {
      const response = await this.drugDiscoveryStore.registerDockingJob(this.experimentId, target.target_id, ligand.ligand_id, target.folding_method!);
      if (response && response.job_id) {
        ligand.jobs.push({ job_id: response.job_id });
      }
    },
    async deleteJob(target: ExtendedTargetMetaData, ligand: ExtendedLigandMetaData, job_id: string, jobIndex: number) {
      const response = await this.drugDiscoveryStore.deleteDockingJob(this.experimentId, target.target_id, ligand.ligand_id, job_id);
      if (response) {
        ligand.jobs.splice(jobIndex, 1); // Remove the job from the list
      }
    },
    async getTargetData(target: ExtendedTargetMetaData) {
      target.loadingTargetData = true;
      try {
        target.data = await this.drugDiscoveryStore.fetchTargetData(this.experimentId, target.target_id);
      } catch (error) {
        console.error('Error fetching target data:', error);
      } finally {
        target.loadingTargetData = false;
      }
    },
    async getLigandDataForTarget(target: ExtendedTargetMetaData, ligand: ExtendedLigandMetaData) {
      ligand.loadingLigandData = true;
      try {
        const data = await this.drugDiscoveryStore.fetchLigandDataForTarget(this.experimentId, target.target_id, ligand.ligand_id);
        ligand.data = {
          sdf_file: data?.ligandSdf,
          smiles: data?.ligandSmiles,
        };
        const jobs = await this.drugDiscoveryStore.getDockingJobsListForTargetLigand(this.experimentId, target.target_id, ligand.ligand_id);
        ligand.jobs = jobs?.jobs_list as JobMetaData[];// Store the jobs in the ligand object
      } catch (error) {
        console.error('Error fetching ligand data:', error);
      } finally {
        ligand.loadingLigandData = false;
      }
    },
    async deleteTarget(targetToDelete: TargetMetaData) {
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: `Loading ligands`
      });
      try {
        await this.drugDiscoveryStore.deleteTargetFromExperiment(this.experimentId, targetToDelete.target_id);
        this.targets = this.targets.filter(target => target.target_id !== targetToDelete.target_id);
      } catch (error) {
        console.error('Error deleting target:', error);
      }
      finally {
        this.$q.loading.hide();
      }
    },
    async deleteLigand(target: ExtendedTargetMetaData, ligandToDelete: ExtendedLigandMetaData) {
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: `Deleting ligand`
      });
      try {
        await this.drugDiscoveryStore.deleteLigandFromTarget(this.experimentId, target.target_id, ligandToDelete.ligand_id);
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
