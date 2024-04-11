<template>
  <q-card v-if="currentJob" flat bordered class="q-pa-sm">
    <div class="row">
      <div class="col-8">
        <q-item-label class="text-h6 q-pa-sm">Current job: {{ currentJob.job_id }}</q-item-label>
        <div class="q-pa-sm">
          <div>Target: {{ currentJob.target_name }}</div>
          <div>Ligand: {{ currentJob.ligand_name }}</div>
          <div>Folding Method: {{ currentJob.folding_method }}</div>
          <div>Docking Method: {{ currentJob.docking_method }}</div>
        </div>
        <q-card-section>
          Prediction steps:
          <div class="q-gutter-md q-pt-sm row items-start">
            <q-card v-for="step in jobSteps" :key="step.label" flat bordered class="my-card">
              <q-card-section :class="getStepClass(currentJob, step.key)">
                <div class="step-content">
                  {{ step.label }}
                  <q-icon v-if="!getServiceHealth(step.key)" name="warning" color="negative" class="text-warning"/>
                  <q-tooltip v-if="!getServiceHealth(step.key)">
                    The {{ step.label }} service is not responding. Check if it's running.
                  </q-tooltip>
                  <template v-if="currentJob[`${step.key}Available`]">
                    <q-icon name="check" class="check-mark" color="green"/>
                  </template>
                  <template v-else-if="currentJob[`${step.key}Running`]">
                    <q-spinner class="spinner" color="primary" size="20px"/>
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
            It appears that predictions are not currently running for this job. Press the "Run" button to continue
            docking.
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
                  :disable="currentJob.isAnyJobRunning || isAnyServiceUnhealthy"
              />
            <q-tooltip v-if="isAnyServiceUnhealthy">
              Fix the services which are not responding to run the docking. Services involved:
            </q-tooltip>
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
        <q-tr :props="props" @click="() => fetchResultDataForRow(props.row)">
          <q-td
              v-for="col in resultsColumns.filter(c => c.name !== 'delete')"
              :key="col.name"
              :props="props"
          >
            {{ props.row[col.field] }}
          </q-td>
          <q-td auto-width>
            <q-btn icon="delete" color="negative" flat @click.stop="deleteDockingJob(props.row)"/>
          </q-td>
        </q-tr>
        <q-tr v-show="props.row.expanded" :props="props">
          <q-td colspan="100%">
            <div class="q-pa-md">
              <template v-if="props.row.docking_method === 'diffdock'">
                <DiffDockResult
                    :experimentId="experimentId"
                    :targetId="props.row.target_id"
                    :ligandId="props.row.ligand_id"
                    :jobId="props.row.job_id"
                />
              </template>
              <template v-else-if="props.row.docking_method === 'umol'">
                <UmolResult
                    :experimentId="experimentId"
                    :targetId="props.row.target_id"
                    :ligandId="props.row.ligand_id"
                    :jobId="props.row.job_id"
                />
              </template>
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
      <q-tr :props="props">
        <q-td v-for="col in queueColumns" :key="col.name" :props="props">
          <template v-if="col.field === 'folding_method'">
            <q-select filled v-model="props.row.folding_method" :options="Object.values(foldingMethods)"
                      @update:modelValue="handleFoldingChange(props.row.job_id, props.row.folding_method)"/>
          </template>
          <template v-else-if="col.field === 'docking_method'">
            <q-select filled v-model="props.row.docking_method" :options="['diffdock', 'umol']"
                      @update:modelValue="handleDockingChange(props.row.job_id, props.row.docking_method)"/>
            <template v-if="props.row.docking_method === 'diffdock'">
              <q-input filled type="number" v-model="props.row.samples_per_complex" label="Samples per Complex"
                       @change="() => updateDiffDockParams(props.row, props.row.samples_per_complex)"/>
            </template>
            <template v-else-if="props.row.docking_method === 'umol'">
              <!-- Placeholder for Pocket IDs Input Field -->
              <div>Pocket IDs input will go here</div>
            </template>
          </template>
          <template v-else-if="col.field === 'progress' && props.row.docking_method === 'diffdock'">
            {{ props.row.progress }}/2
          </template>
          <template v-else-if="col.field === 'progress' && props.row.docking_method === 'umol'">
            {{ props.row.progress }}/4
          </template>
          <template v-else-if="col.field === 'delete'">
            <q-btn icon="delete" color="negative" flat @click.stop="deleteDockingJob(props.row)"/>
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
              :disabled="isAnyServiceUnhealthy || isAnyJobRunning"
          />
          <q-tooltip v-if="isAnyServiceUnhealthy">
            Fix the services which are not responding to run the docking. Services involved
          </q-tooltip>
          <q-btn icon="expand_more" @click.stop="props.expand = !props.expand"/>
        </q-td>
      </q-tr>
      <q-tr v-if="props.expand" :props="props">
        <q-td colspan="100%">
          <!-- Expanded row content -->
          <div v-if="props.row.pocketIds">
            <div><strong>Pocket IDs:</strong> {{ props.row.pocketIds.join(', ') }}</div>
          </div>
          <div v-else>
            Binding pocket was not set manually, will be predicted
          </div>
        </q-td>
      </q-tr>
    </template>
  </q-table>
</template>

<script lang="ts">
import {useDrugDiscoveryStore} from 'src/features/drug_discovery/storage';
import {useRoute} from "vue-router";
import DiffDockResult from "src/features/drug_discovery/components/results/DiffDockResult.vue";
import UmolResult from "src/features/drug_discovery/components/results/UmolResult.vue";
import {defineComponent} from "vue";
import {JobMetaData} from "src/api/client";

interface Job {
  job_id: string;
  target_name: string | undefined;
  ligand_name: string | undefined;
  target_id: string;
  ligand_id: string;
  folding_method: string;
  docking_method: string;
  progress?: number;
  isRunning?: boolean;
  msaAvailable?: boolean;
  pocketAvailable?: boolean;
  foldingAvailable?: boolean;
  dockingAvailable?: boolean;
  msaRunning?: boolean;
  pocketRunning?: boolean;
  foldingRunning?: boolean;
  dockingRunning?: boolean;
  pocketIds: number[] | null;
  isAnyJobRunning: boolean;
  samples_per_complex?: number;

  [key: string]: any; // Add this line to allow dynamic keys
}

interface ResultData {
  predicted_pdb?: string;
  predicted_sdf?: string;
  predicted_pdb_file?: File;
  predicted_sdf_file?: File;
}

interface DockingResult {
  target_name: string | undefined;
  ligand_name: string | undefined;
  folding_method: string;
  docking_method: string;
  resultData?: ResultData | null;
  experiment_id: string;
  job_id: string;
  target_id: string;
  ligand_id: string;
  expanded: boolean;
}

export default defineComponent({
  name: 'RunDocking',
  components: {DiffDockResult, UmolResult},
  data() {
    return {
      foldingMethods: {
        esmfoldLight: 'esmfold_light',
        esmfold: 'esmfold',
        rosettafold: 'rosettafold'
      },
      dockingResults: [] as DockingResult[],
      jobsInQueue: [] as Job[],
      resultsColumns: [
        {name: 'job_id', label: 'Job ID', field: 'job_id', align: 'left'},
        {name: 'target_name', label: 'Target Name', field: 'target_name', align: 'left'},
        {name: 'ligand_name', label: 'Ligand Name', field: 'ligand_name', align: 'left'},
        {name: 'folding_method', label: 'Folding Method', field: 'folding_method', align: 'center'},
        {name: 'docking_method', label: 'Docking Method', field: 'docking_method', align: 'center'},
        {name: 'delete', label: 'Delete', field: 'delete', align: 'center'},
      ],
      queueColumns: [
        {name: 'job_id', label: 'Job ID', field: 'job_id', align: 'left'},
        {name: 'target_name', label: 'Target Name', field: 'target_name', align: 'left'},
        {name: 'ligand_name', label: 'Ligand Name', field: 'ligand_name', align: 'left'},
        {name: 'folding_method', label: 'Folding Method', field: 'folding_method', align: 'center'},
        {name: 'docking_method', label: 'Docking Method', field: 'docking_method', align: 'center'},
        {name: 'progress', label: 'Progress', field: 'progress', align: 'left', sortable: true},
        {name: 'delete', label: 'Delete', field: 'delete', align: 'center'},
      ],
      currentJob: null as Job | null,
      msaServiceHealthy: false,
      p2RankServiceHealthy: false,
      foldingServiceHealthy: false,
      dockingServiceHealthy: false,
      polling: null as number | null,
    };
  },
  methods: {
    async updateServiceHealth() {
      const store = useDrugDiscoveryStore();
      const msaServiceHealthyResponse = await store.checkMsaServiceHealth();
      const p2RankServiceHealthyResponse = await store.checkP2RankServiceHealth();

      this.msaServiceHealthy = msaServiceHealthyResponse?.is_healthy ?? false;
      this.p2RankServiceHealthy = p2RankServiceHealthyResponse?.is_healthy ?? false;

      if (this.currentJob?.docking_method === 'diffdock') {
        const diffdockServiceHealthyResponse = await store.checkDiffDockServiceHealth();
        this.dockingServiceHealthy = diffdockServiceHealthyResponse?.is_healthy ?? false;
      } else if (this.currentJob?.docking_method === 'umol') {
        const umolServiceHealthyResponse = await store.checkUmolServiceHealth();
        this.dockingServiceHealthy = umolServiceHealthyResponse?.is_healthy ?? false;
      }
      // Check folding service health based on the folding method of the current job
      if (this.currentJob?.folding_method === this.foldingMethods.esmfold) {
        const healthy = await store.checkEsmFoldServiceHealth();
        this.foldingServiceHealthy = healthy?.is_healthy ?? false;
      }
      if (this.currentJob?.folding_method === this.foldingMethods.esmfoldLight) {
        const healthy = await store.checkEsmFoldLightServiceHealth();
        this.foldingServiceHealthy = healthy?.is_healthy ?? false;
      }
      if (this.currentJob?.folding_method === this.foldingMethods.rosettafold) {
        const healthy = await store.checkRosettaFoldServiceHealth();
        this.foldingServiceHealthy = healthy?.is_healthy ?? false;
      }
    },
    getServiceHealth(key: string) {
      switch (key) {
        case 'msa':
          return this.msaServiceHealthy;
        case 'pocket':
          return this.p2RankServiceHealthy;
        case 'folding':
          return this.foldingServiceHealthy;
        case 'docking':
          return this.dockingServiceHealthy;
        default:
          return true; // Default to healthy if the service key is not recognized
      }
    },
    async fetchDockingResults() {
      // experimentId is definitely not null here so can use "!"
      const store = useDrugDiscoveryStore();
      const results = await store.getAllDockingResultsList(this.experimentId!);

      for (const result of results?.results_list ?? []) {
        const targetMeta = await store.fetchTargetMetaData(this.experimentId!, result.target_id);
        const ligandMeta = await store.fetchLigandMetaDataForTarget(this.experimentId!, result.target_id, result.ligand_id);

        const dockingResult = {
          ...result,
          target_name: targetMeta?.target_name,
          ligand_name: ligandMeta?.ligand_name,
          resultData: null, // Initialize resultData as null
          expanded: false,
        };
        this.dockingResults.push(dockingResult);
      }
    },

    async populateDiffDockParamsForJobs() {
      for (let job of this.jobsInQueue) {
        if (job.docking_method === 'diffdock') {
          await this.loadDiffDockParams(job);
        }
      }
    },

    async loadDiffDockParams(job: Job) {
      const store = useDrugDiscoveryStore();
      const params = await store.getDiffDockParams(this.experimentId!, job.target_id, job.ligand_id, job.job_id);
      if (params && params.samples_per_complex !== undefined) {
        job.samples_per_complex = params.samples_per_complex;
      }
    },
    async updateDiffDockParams(job: JobMetaData, samplesPerComplex: number) {
      const store = useDrugDiscoveryStore();
      await store.updateDiffDockParams(this.experimentId!, job.target_id, job.ligand_id, job.job_id, samplesPerComplex);
    },

    async deleteDockingJob(row: DockingResult) {
      const store = useDrugDiscoveryStore();
      // Call the API to delete the docking job
      await store.deleteDockingJob(this.experimentId!, row.target_id, row.ligand_id, row.job_id);

      // Remove the job from the dockingResults array
      this.dockingResults = this.dockingResults.filter(job => job.job_id !== row.job_id);

      // Remove the job from the jobsInQueue array
      this.jobsInQueue = this.jobsInQueue.filter(job => job.job_id !== row.job_id);

      // You might want to check if the currentJob needs to be updated
      if (this.currentJob && this.currentJob.job_id === row.job_id) {
        this.currentJob = null; // or set it to another job as per your logic
      }
    },

    async fetchResultDataForRow(row: DockingResult) {
      row.expanded = !row.expanded;
    },

    async fetchJobsInQueue() {
      const store = useDrugDiscoveryStore();
      const jobs = await store.getAllDockingJobsList(this.experimentId!);

      for (const job of jobs?.jobs_list ?? []) {
        const targetMeta = await store.fetchTargetMetaData(this.experimentId!, job.target_id);
        const ligandMeta = await store.fetchLigandMetaDataForTarget(this.experimentId!, job.target_id, job.ligand_id);
        const pocketIdsResponse = await store.getJobPocketIds(this.experimentId!, job.target_id, job.ligand_id, job.job_id);

        let pocketIds = null;
        if (pocketIdsResponse && pocketIdsResponse.pocket_ids) {
          pocketIds = pocketIdsResponse.pocket_ids;
        }
        this.jobsInQueue.push({
          ...job,
          target_name: targetMeta?.target_name,
          ligand_name: ligandMeta?.ligand_name,
          progress: await this.calculateProgress(job.target_id, job.ligand_id, job.job_id, job.docking_method, job.folding_method),
          pocketIds: pocketIds || null,
          isAnyJobRunning: false
        });
      }
    },

    async calculateProgress(target_id: string, ligand_id: string, job_id: string, docking_method: string, folding_method: string) {
      const store = useDrugDiscoveryStore();
      let progress = 0;

      const foldingAvailable = await store.checkFoldingDataAvailable(this.experimentId!, target_id, folding_method);
      if (foldingAvailable?.is_available) progress++;

      const dockingAvailable = await store.checkDockingResultAvailable(this.experimentId!, target_id, ligand_id, job_id);
      if (dockingAvailable?.result_available) progress++;

      // For diffdock, only consider folding and docking results
      if (docking_method === 'diffdock') {
        return progress;
      }

      // For other docking methods, consider all steps (including MSA and pocket prediction)
      const msaAvailable = await store.checkMsaDataAvailable(this.experimentId!, target_id);
      if (msaAvailable?.is_available) progress++;

      const pocketAvailable = await store.checkPocketDataAvailable(this.experimentId!, target_id, ligand_id, job_id);
      if (pocketAvailable?.is_available) progress++;

      return progress;
    },

    runCurrentJob() {
      if (!this.currentJob?.isAnyJobRunning) {
        if (this.currentJob != null)
          this.runDockingJob(this.currentJob);
      }
    },

    async runDockingJob(row: Job) {
      if (this.isAnyJobRunning) return;

      // Mark the selected job as running
      this.jobsInQueue = this.jobsInQueue.map(job => ({
        ...job,
        isRunning: job.job_id === row.job_id
      }));

      // Call the store method to start the job based on the docking method
      const store = useDrugDiscoveryStore();
      if (row.docking_method === 'diffdock') {
        await store.runDiffDockDockingJob(this.experimentId!, row.target_id, row.ligand_id, row.job_id);
      } else if (row.docking_method === 'umol') {
        await store.runUmolDockingJob(this.experimentId!, row.target_id, row.ligand_id, row.job_id);
      }

    },

    async handleFoldingChange(job_id: string, folding_method: string) {
      const store = useDrugDiscoveryStore();
      const job = this.jobsInQueue.find(j => j.job_id === job_id);
      if (job) {
        await store.updateDockingParams(this.experimentId!, job.target_id, job.ligand_id, job.job_id, folding_method, job.docking_method);
        job.folding_method = folding_method;
        job.progress = await this.calculateProgress(job.target_id, job.ligand_id, job.job_id, job.docking_method, job.folding_method);
      }
    },
    async handleDockingChange(job_id: string, docking_method: string) {
      const store = useDrugDiscoveryStore();
      const job = this.jobsInQueue.find(j => j.job_id === job_id);
      if (job) {
        await store.updateDockingParams(this.experimentId!, job.target_id, job.ligand_id, job.job_id, job.folding_method, docking_method);
        job.docking_method = docking_method;
        job.progress = await this.calculateProgress(job.target_id, job.ligand_id, job.job_id, job.docking_method, job.folding_method);
      }
    },

    async updateCurrentJob() {
      const store = useDrugDiscoveryStore();
      await this.updateServiceHealth(); // Update the health status of each service
      const runningJob = this.jobsInQueue.find(job => job.isRunning);
      this.currentJob = runningJob || this.jobsInQueue[0] || null;

      if (!this.currentJob) return;

      const msaAvailableResp = await store.checkMsaDataAvailable(this.experimentId!, this.currentJob.target_id);
      this.currentJob.msaAvailable = msaAvailableResp?.is_available ?? false;
      if (!this.currentJob.msaAvailable && this.msaServiceHealthy) {
        const msaRunningResp = await store.checkMsaJobIsRunning(this.currentJob.job_id);
        this.currentJob.msaRunning = msaRunningResp?.is_running ?? false;
      }

      const pocketAvailableResp = await store.checkPocketDataAvailable(this.experimentId!, this.currentJob.target_id, this.currentJob.ligand_id, this.currentJob.job_id);
      this.currentJob.pocketAvailable = pocketAvailableResp?.is_available ?? false;
      if (!this.currentJob.pocketAvailable && this.p2RankServiceHealthy) {
        const p2RankRunningResp = await store.checkP2RankJobIsRunning(this.currentJob.job_id);
        this.currentJob.pocketRunning = p2RankRunningResp?.is_running ?? false;
      }

      const foldingAvailableResp = await store.checkFoldingDataAvailable(this.experimentId!, this.currentJob.target_id, this.currentJob.folding_method);
      this.currentJob.foldingAvailable = foldingAvailableResp?.is_available ?? false;
      if (!this.currentJob.foldingAvailable && this.foldingServiceHealthy) {
        let foldingRunningResp;
        if (this.currentJob.folding_method === this.foldingMethods.esmfold) {
          foldingRunningResp = await store.checkEsmFoldJobIsRunning(this.currentJob.job_id);
        }
        if (this.currentJob.folding_method === this.foldingMethods.esmfoldLight) {
          foldingRunningResp = await store.checkEsmFoldLightJobIsRunning(this.currentJob.job_id);
        }
        if (this.currentJob.folding_method === this.foldingMethods.rosettafold) {
          foldingRunningResp = await store.checkRosettafoldJobIsRunning(this.currentJob.job_id);
        }
        this.currentJob.foldingRunning = foldingRunningResp?.is_running ?? false;
      }

      const dockingResultAvailableResp = await store.checkDockingResultAvailable(this.experimentId!, this.currentJob.target_id, this.currentJob.ligand_id, this.currentJob.job_id);
      this.currentJob.dockingAvailable = dockingResultAvailableResp?.result_available ?? false;
      if (!this.currentJob.dockingAvailable && this.dockingServiceHealthy) {
        let dockingRunningResp;
        if (this.currentJob.docking_method === 'diffdock') {
          dockingRunningResp = await store.checkDiffDockJobIsRunning(this.currentJob.job_id);
        } else if (this.currentJob.docking_method === 'umol') {
          dockingRunningResp = await store.checkUmolJobIsRunning(this.currentJob.job_id);
        }
        this.currentJob.dockingRunning = dockingRunningResp?.is_running ?? false;
      }

      // Update the overall running status based on individual service running status
      if (this.currentJob) {
        this.currentJob.isAnyJobRunning = (this.currentJob.msaRunning || this.currentJob.pocketRunning || this.currentJob.foldingRunning || this.currentJob.dockingRunning) ?? false;
      }

      let allDataAvailable;

      if (this.currentJob.docking_method === 'diffdock') {
        allDataAvailable = this.currentJob.foldingAvailable && this.currentJob.dockingAvailable;
      } else if (this.currentJob.docking_method === 'umol') {
        allDataAvailable = this.currentJob.msaAvailable && this.currentJob.pocketAvailable && this.currentJob.foldingAvailable && this.currentJob.dockingAvailable;
      }

      // If all data is available, refresh the current job and results
      if (allDataAvailable) {
        this.dockingResults = [];
        this.jobsInQueue = [];

        // Refresh the docking results and the job queue
        await this.fetchDockingResults();
        await this.fetchJobsInQueue();
      }
    },
    getStepClass(job: Job | null, step: string) {
      if (!job) return 'bg-grey-3 text-black q-ma-md';

      if (job[`${step}Available`]) {
        return 'bg-green-3 text-black q-ma-md';
      } else if (job[`${step}Running`]) {
        return 'bg-grey-3 text-black q-ma-md';
      } else {
        return 'bg-grey-3 text-black q-ma-md';
      }
    },
  },
  computed: {
    experimentId() {
      const route = useRoute();
      return route.params.experimentId as string;
    },
    isAnyServiceUnhealthy() {
      if (this.currentJob?.docking_method === 'diffdock') {
        return  !this.foldingServiceHealthy || !this.dockingServiceHealthy;
      } else if (this.currentJob?.docking_method === 'umol') {
        return !this.msaServiceHealthy || !this.p2RankServiceHealthy || !this.foldingServiceHealthy || !this.dockingServiceHealthy;
      }
      console.error('No docking method selected');
      return false;
      // Computed property to check if any service is unhealthy
    },
    isAnyJobRunning() {
      return this.jobsInQueue.some(job => job.isRunning);
    },
    jobSteps() {
      if (this.currentJob?.docking_method === 'diffdock') {
        return [
          {key: 'folding', label: 'Step 1. Folding'},
          {key: 'docking', label: 'Step 2. Docking'}
        ];
      } else {
        return [
          {key: 'msa', label: 'Step 1. MSA'},
          {key: 'pocket', label: 'Step 2. Binding Pocket'},
          {key: 'folding', label: 'Step 3. Folding'},
          {key: 'docking', label: 'Step 4. Docking'},
        ];
      }
    }
  },
  async mounted() {
    this.polling = window.setInterval(this.updateCurrentJob, 1000);
    await this.fetchDockingResults();
    await this.fetchJobsInQueue();
    await this.populateDiffDockParamsForJobs();
  },
  unmounted() {
    if (this.polling !== null) {
      clearInterval(this.polling); // TypeScript knows this.polling is not null here
      this.polling = null;
    }
  }
})
</script>
