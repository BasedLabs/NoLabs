<template>
  <div>
    <q-card v-if="currentJob" flat bordered class="q-pa-sm">
      <q-item-label v-if="currentJob" class="text-h6 q-pa-sm"> Current job running: {{ currentJob.job_id }} </q-item-label>
      <q-item-label v-if="currentJob" class="text-h8 q-pa-sm"> TargetId: {{ currentJob.target_id }} </q-item-label>
      <q-item-label v-if="currentJob" class="text-h8 q-pa-sm"> LigandId: {{ currentJob.ligand_id }} </q-item-label>
      <q-card-section>
        <q-item-label v-if="currentJob" class="text-h8 q-pa-sm"> Prediction steps: </q-item-label>
        <div class="q-gutter-md q-pt-sm row items-start">
          <q-card v-for="step in steps" :key="step.label" flat bordered class="my-card">
            <q-card-section :class="getStepClass(currentJob, step.key)">
              <div class="step-content">
                {{ step.label }}
                <template v-if="currentJob && currentJob[`${step.key}Running`]">
                  <q-spinner class="spinner" color="primary" size="20px" />
                </template>
                <template v-else-if="currentJob && currentJob[`${step.key}Available`]">
                  <q-icon name="check" class="check-mark" color="info" />
                </template>
              </div>
            </q-card-section>
          </q-card>
        </div>
      </q-card-section>
    </q-card>

    <q-card flat bordered class="q-pa-sm q-mt-md q-mb-md">
      <q-table
        title="Docking Results"
        :rows="dockingResults"
        :columns="resultsColumns"
        row-key="job_id"
        class="q-pa-md"
      />
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

      <template v-slot:body-cell-action="props">
        <q-td :props="props">
          <q-btn
            label="Run"
            color="primary"
            @click="runDockingJob(props.row)"
            :disabled="props.row.isRunning"
          />
        </q-td>
      </template>
    </q-table>
  </div>
</template>

<script>
import { useDrugDiscoveryStore } from 'src/features/drug_discovery/storage';
import {useRoute} from "vue-router";

export default {
  name: 'RunDocking',
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
        { name: 'action', label: 'Action', field: 'action', align: 'left' },
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
        const isDataAvailableResponse = await store.checkDockingResultAvailable(this.experimentId,
          result.target_id,
          result.ligand_id,
          result.job_id);
        if (isDataAvailableResponse.result_available) {
          const targetMeta = await store.fetchTargetMetaData(this.experimentId, result.target_id);
          const ligandMeta = await store.fetchLigandMetaData(this.experimentId, result.target_id, result.ligand_id);
          this.dockingResults.push({
            ...result,
            target_name: targetMeta.target_name,
            ligand_name: ligandMeta.ligand_name,
          });
        }
      }
    },

    async fetchJobsInQueue() {
      const store = useDrugDiscoveryStore();
      const results = await store.getAllDockingResultsList(this.experimentId);

      for (const result of results.results_list) {
        const targetMeta = await store.fetchTargetMetaData(this.experimentId, result.target_id);
        const ligandMeta = await store.fetchLigandMetaData(this.experimentId, result.target_id, result.ligand_id);
        const isRunning = await this.checkIfAnyJobIsRunning(result.job_id);

        this.jobsInQueue.push({
          ...result,
          target_name: targetMeta.target_name,
          ligand_name: ligandMeta.ligand_name,
          progress: await this.calculateProgress(result),
          isRunning,
        });
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

    async checkIfAnyJobIsRunning(jobId) {
      const store = useDrugDiscoveryStore();
      const msaRunning = await store.checkMsaJobIsRunning(jobId);
      const foldingRunning = await store.checkFoldingJobIsRunning(jobId);
      const p2RankRunning = await store.checkP2RankJobIsRunning(jobId);
      const umolRunning = await store.checkUmolJobIsRunning(jobId);

      return msaRunning.is_running || foldingRunning.is_running || p2RankRunning.is_running || umolRunning.is_running;
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
      }
    },

    getStepClass(job, step) {
      if (!job) return 'bg-grey-3 text-black q-ma-md';

      if (job[`${step}Running`]) {
        return 'bg-primary text-white q-ma-md';
      } else if (job[`${step}Available`]) {
        return 'bg-green-3 text-black q-ma-md';
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
