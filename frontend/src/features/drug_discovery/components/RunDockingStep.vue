<template>
    <q-card v-if="currentJob" flat bordered class="q-pa-sm">
      <div class="row">
        <div class="col-8">
          <q-item-label class="text-h6 q-pa-sm">Current job: {{ currentJob.job_id }}</q-item-label>
          <div class="q-pa-sm">
            <div>Target: {{ currentJob.target_name }}</div>
            <div>Ligand: {{ currentJob.ligand_name }}</div>
          </div>
          <q-card-section>
            Prediction steps:
            <div class="q-gutter-md q-pt-sm row items-start">
              <q-card v-for="step in steps" :key="step.label" flat bordered class="my-card">
                <q-card-section :class="getStepClass(currentJob, step.key)">
                  <div class="step-content">
                    {{ step.label }}
                    <template v-if="currentJob[`${step.key}Available`]">
                      <q-icon name="check" class="check-mark" color="green" />
                    </template>
                    <template v-else-if="currentJob[`${step.key}Running`]">
                      <q-spinner class="spinner" color="primary" size="20px" />
                    </template>
                  </div>
                </q-card-section>
              </q-card>
            </div>
          </q-card-section>
        </div>
        <div class="col-4">
          <div class="q-pa-md col items-center">
            <div v-if="!currentJob.isAnyJobRunning" class="text-warning q-pa-md">
              It appears that predictions are not currently running for this job. Press the "Run" button to continue docking.
            </div>
            <div class="q-pa-md" v-else>
              Job is running...
            </div>
            <div v-if="!currentJob.isAnyJobRunning" class="q-pa-md">
              <q-btn
                label="Run current job"
                color="info"
                @click="runCurrentJob"
                :loading="currentJob.isAnyJobRunning"
                :disable="currentJob.isAnyJobRunning"
              />
            </div>
            <div class="q-pa-md" v-else>
              <q-spinner size="xl"></q-spinner>
            </div>
          </div>
        </div>
      </div>
    </q-card>

  <q-card flat bordered class="q-pa-sm q-mt-md q-mb-md">
    <q-table
      title="Docking Results"
      :rows="dockingResults"
      :columns="resultsColumns"
      row-key="job_id"
      class="q-pa-md bg-black"
    >
      <template v-slot:body="props">
        <q-tr :props="props" @click.stop="fetchResultDataForRow(props.row)">
          <q-td v-for="col in resultsColumns" :key="col.name" :props="props">
            {{ props.row[col.field] }}
          </q-td>
        </q-tr>
        <q-tr v-if="props.row.expanded" :props="props">
          <q-td colspan="100%">
            <div class="q-pa-md">
              <div v-if="props.row.resultData">
                <!-- Display docking job result data -->
                <PdbViewer v-if="props.row.resultData.predicted_sdf_file" :pdb-file="props.row.resultData.predicted_pdb_file" :sdf-file="props.row.resultData.predicted_sdf_file"/>
                <div><strong>PLDDT Array:</strong> {{ props.row.resultData.plddt_array.join(', ') }}</div>
              </div>
              <div v-else>
                <q-spinner color="info" label="Loading..." />
              </div>
            </div>
          </q-td>
        </q-tr>
      </template>
    </q-table>
  </q-card>


  <q-table
      title="Jobs in Queue"
      :rows="jobsInQueue"
      :columns="queueColumns"
      row-key="job_id"
    >
      <template v-slot:body-cell-progress="props">
        <q-td :props="props">
          {{ props.row.progress }}/4
        </q-td>
      </template>

      <template v-slot:body="props">
        <q-tr :props="props" @click.stop="props.expand = !props.expand">
          <q-td v-for="col in queueColumns" :key="col.name" :props="props">
            <template v-if="col.field === 'progress'">
              {{ props.row.progress }}/4
            </template>
            <template v-else>
              {{ props.row[col.field] }}
            </template>
          </q-td>
          <q-td auto-width>
            <q-btn
              label="Run"
              color="info"
              @click.stop="runDockingJob(props.row)"
              :disabled="props.row.isRunning"
            />
            <q-btn icon="expand_more" @click.stop="props.expand = !props.expand" />
          </q-td>
        </q-tr>
        <q-tr v-if="props.expand" :props="props">
          <q-td colspan="100%">
            <!-- Expanded row content -->
            Job input parameters for {{ props.row.job_id }}
          </q-td>
        </q-tr>
      </template>
    </q-table>
</template>

<script>
import { useDrugDiscoveryStore } from 'src/features/drug_discovery/storage';
import {useRoute} from "vue-router";
import PdbViewer from "src/components/PdbViewer.vue";

export default {
  name: 'RunDocking',
  components: {PdbViewer},
  data() {
    return {
      dockingResults: [],
      jobsInQueue: [],
      resultsColumns: [
        { name: 'job_id', label: 'Job ID', field: 'job_id', align: 'left' },
        { name: 'target_name', label: 'Target Name', field: 'target_name', align: 'left' },
        { name: 'ligand_name', label: 'Ligand Name', field: 'ligand_name', align: 'left' },
      ],
      queueColumns: [
        { name: 'job_id', label: 'Job ID', field: 'job_id', align: 'left' },
        { name: 'target_name', label: 'Target Name', field: 'target_name', align: 'left' },
        { name: 'ligand_name', label: 'Ligand Name', field: 'ligand_name', align: 'left' },
        { name: 'progress', label: 'Progress', field: 'progress', align: 'left', sortable: true },
      ],
      currentJob: null,
      steps: [
        { key: 'msa', label: 'Step 1. MSA' },
        { key: 'pocket', label: 'Step 2. Binding Pocket' },
        { key: 'folding', label: 'Step 3. Folding' },
        { key: 'docking', label: 'Step 4. Docking' },
      ],
    };
  },
  methods: {
    async fetchDockingResults() {
      const store = useDrugDiscoveryStore();
      const results = await store.getAllDockingResultsList(this.experimentId);

      for (const result of results.results_list) {
        const targetMeta = await store.fetchTargetMetaData(this.experimentId, result.target_id);
        const ligandMeta = await store.fetchLigandMetaData(this.experimentId, result.target_id, result.ligand_id);

        const dockingResult = {
          ...result,
          target_name: targetMeta.target_name,
          ligand_name: ligandMeta.ligand_name,
          resultData: null, // Initialize resultData as null
        };

        this.dockingResults.push(dockingResult);
      }
    },

    async fetchResultDataForRow(row) {
      if (!row.expanded && !row.resultData) {
        const store = useDrugDiscoveryStore();
        row.resultData = await store.getDockingJobResultData(this.experimentId, row.target_id, row.ligand_id, row.job_id);
        row.resultData.predicted_pdb_file = new File([new Blob([row.resultData.predicted_pdb])], "protein.pdb");
        row.resultData.predicted_sdf_file = new File([new Blob([row.resultData.predicted_sdf])], "ligand.sdf");
      }
      row.expanded = !row.expanded; // Toggle row expansion
    },

    async fetchJobsInQueue() {
      const store = useDrugDiscoveryStore();
      const results = await store.getAllDockingResultsList(this.experimentId);

      for (const result of results.results_list) {
        const targetMeta = await store.fetchTargetMetaData(this.experimentId, result.target_id);
        const ligandMeta = await store.fetchLigandMetaData(this.experimentId, result.target_id, result.ligand_id);
        const isDataAvailableResponse = await store.checkDockingResultAvailable(this.experimentId,
          result.target_id,
          result.ligand_id,
          result.job_id);
        if (!isDataAvailableResponse.result_available) {
          this.jobsInQueue.push({
            ...result,
            target_name: targetMeta.target_name,
            ligand_name: ligandMeta.ligand_name,
            progress: await this.calculateProgress(result)
          });
        }
      }
    },

    async calculateProgress(result) {
      const store = useDrugDiscoveryStore();
      let progress = 0;

      const msaAvailable = await store.checkMsaDataAvailable(this.experimentId, result.target_id);
      if (msaAvailable.is_available) progress++;

      const foldingAvailable = await store.checkFoldingDataAvailable(this.experimentId, result.target_id);
      if (foldingAvailable.is_available) progress++;

      const pocketAvailable = await store.checkPocketDataAvailable(this.experimentId, result.target_id, result.ligand_id, result.job_id);
      if (pocketAvailable.is_available) progress++;

      const dockingAvailable = await store.checkDockingResultAvailable(this.experimentId, result.target_id, result.ligand_id, result.job_id);
      if (dockingAvailable.result_available) progress++;

      return progress;
    },

    runCurrentJob() {
      if (!this.currentJob.isAnyJobRunning) {
        this.runDockingJob(this.currentJob);
      }
    },

    async runDockingJob(row) {
      const store = useDrugDiscoveryStore();
      await store.runDockingJob(this.experimentId, row.target_id, row.ligand_id, row.job_id);
      // You might want to refresh the job list or perform some other action upon starting the job.
    },

    async updateCurrentJob() {
      const store = useDrugDiscoveryStore();
      this.currentJob = this.jobsInQueue.length ? this.jobsInQueue[0] : null;

      if (this.currentJob) {
        // Check if data is available for each step
        const msaAvailableResp = await store.checkMsaDataAvailable(this.experimentId, this.currentJob.target_id);
        this.currentJob.msaAvailable = msaAvailableResp.is_available;

        const foldingAvailableResp = await store.checkFoldingDataAvailable(this.experimentId, this.currentJob.target_id);
        this.currentJob.foldingAvailable = foldingAvailableResp.is_available;

        const pocketAvailableResp = await store.checkPocketDataAvailable(this.experimentId, this.currentJob.target_id, this.currentJob.ligand_id, this.currentJob.job_id);
        this.currentJob.pocketAvailable = pocketAvailableResp.is_available;

        const dockingResultAvailableResp = await store.checkDockingResultAvailable(this.experimentId, this.currentJob.target_id, this.currentJob.ligand_id, this.currentJob.job_id);
        this.currentJob.dockingAvailable = dockingResultAvailableResp.result_available;

        // Update the running status for each step
        const msaRunningResp = await store.checkMsaJobIsRunning(this.currentJob.job_id);
        this.currentJob.msaRunning = msaRunningResp.is_running;

        const foldingRunningResp = await store.checkFoldingJobIsRunning(this.currentJob.job_id);
        this.currentJob.foldingRunning = foldingRunningResp.is_running;

        const p2RankRunningResp = await store.checkP2RankJobIsRunning(this.currentJob.job_id);
        this.currentJob.pocketRunning = p2RankRunningResp.is_running;

        const umolRunningResp = await store.checkUmolJobIsRunning(this.currentJob.job_id);
        this.currentJob.dockingRunning = umolRunningResp.is_running;

        this.currentJob.isAnyJobRunning = this.currentJob.dockingRunning || this.currentJob.foldingRunning.is_running || this.currentJob.pocketRunning.is_running || this.currentJob.dockingRunning.is_running;
      }
    },

    getStepClass(job, step) {
      if (!job) return 'bg-grey-3 text-black q-ma-md';

      if (job[`${step}Available`]) {
        return 'bg-green-3 text-black q-ma-md';
      }
      else if (job[`${step}Running`]) {
        return 'bg-grey-3 text-black q-ma-md';
      } else {
        return 'bg-grey-3 text-black q-ma-md';
      }
    },
  },
  mounted() {
    const route = useRoute();
    this.experimentId = route.params.experimentId;
    this.fetchDockingResults();
    this.fetchJobsInQueue();
    this.polling = setInterval(this.updateCurrentJob, 1000);
  },
  beforeUnmount() {
    clearInterval(this.polling); // Clear polling interval when component is destroyed
  },
};
</script>
