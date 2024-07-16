<template>
  <q-card-section class="job-section q-pa-sm">
    <div class="text-white q-pa-sm text-caption">Jobs queue</div>
    <q-scroll-area v-if="jobs.length > 0" visible :thumbStyle="thumbStyle" :barStyle="barStyle" style="height: 30vh">
      <draggable class="q-pa-sm" v-model="jobs" handle=".drag-handle" @end="updateJobOrder" item-key="job_id">
        <template #item="{ element }">
          <q-item class="bg-grey-3 text-black">
            <q-item-section>
              <q-btn @click="copyContent(element.job_id)" icon="content_copy" label="ID"></q-btn>
            </q-item-section>
            <q-item-section class="col">
              <q-item-label class="text-truncate">{{ element.job_name }}</q-item-label>
              <q-tooltip>{{ element.job_name }}</q-tooltip>
            </q-item-section>
            <q-item-section>
              <div v-if="element.executionStatus === null">No execution status</div>
              <div v-else>
                <q-spinner v-if="element.executionStatus.running" color="primary" size="20px" />
                {{ element.executionStatus.running ? 'Running...' : 'Not running' }}
              </div>
            </q-item-section>
            <q-item-section v-if="jobErrors[element.job_id]">
              <q-item-label class="text-red">{{ jobErrors[element.job_id] }}</q-item-label>
            </q-item-section>
            <q-item-section>
              <q-btn to="" @click="openJob(element)" label="View" dense />
            </q-item-section>
            <q-item-section>
              <q-btn @click="deleteJob(element)" color="negative" label="Delete" dense />
            </q-item-section>
          </q-item>
        </template>
      </draggable>
    </q-scroll-area>
  </q-card-section>

  <q-card-section class="result-section q-pa-sm">
    <div class="text-white q-pa-sm text-caption">Completed jobs</div>
    <q-scroll-area v-if="results.length > 0" visible :thumbStyle="thumbStyle" :barStyle="barStyle" style="height: 30vh">
      <draggable class="q-pa-sm" v-model="results" handle=".drag-handle" item-key="job_id">
        <template #item="{ element }">
          <q-item class="bg-grey-3 text-black">
            <q-item-section>
              <q-btn @click="copyContent(element.job_id)" icon="content_copy" label="ID"></q-btn>
            </q-item-section>
            <q-item-section class="col">
              <q-item-label class="text-truncate">{{ element.job_name }}</q-item-label>
              <q-tooltip>{{ element.job_name }}</q-tooltip>
            </q-item-section>
            <q-item-section v-if="jobErrors[element.job_id]">
              <q-item-label class="text-red">{{ jobErrors[element.job_id] }}</q-item-label>
            </q-item-section>
            <q-item-section>
              <q-btn @click="openJob(element)" label="View" dense />
            </q-item-section>
            <q-item-section>
              <q-btn @click="deleteJob(element)" color="negative" label="Delete" dense />
            </q-item-section>
          </q-item>
        </template>
      </draggable>
    </q-scroll-area>
  </q-card-section>

  <q-card-section class="exception-section q-pa-sm" v-if="lastExceptions.length">
    <div class="text-white q-pa-sm">Last Exceptions</div>
    <q-item class="text-black q-pa-sm q-mb-sm q-border-radius-md">
      <q-item-section>
        <q-btn color="red" label="Show last exceptions" class="q-pm-md"
          @click="showLastExceptionsModal = !showLastExceptionsModal">
          <q-icon left name="warning" />
        </q-btn>
      </q-item-section>
      <ComponentExceptionsModal v-model:visible="showLastExceptionsModal" :exceptions="lastExceptions" />
    </q-item>
  </q-card-section>

  <q-dialog v-model="showJobModal" full-width>
    <q-card style="max-width: 90vw;">
      <component :is="selectedJobComponent" :job-id="selectedJobId" />
      <q-card-actions>
        <q-btn flat label="Close" v-close-popup />
      </q-card-actions>
    </q-card>
  </q-dialog>
</template>

<script lang="ts">
import { defineComponent, watch } from 'vue';
import { useWorkflowStore, Node } from 'src/features/workflow/components/storage';
import {
  getJobMetaData
} from 'src/features/workflow/refinedApi';
import {
  GetJobMetadataResponse,
  nolabs__application__use_cases__folding__api_models__GetJobStatusResponse,
} from 'src/refinedApi/client';
import EsmFoldJob from '../jobs/EsmFoldJob.vue';
import P2RankJob from '../jobs/P2RankJob.vue';
import DiffDockJob from '../jobs/DiffDockJob.vue';
import draggable from 'vuedraggable';
import { Notify } from "quasar";
import componentApi from "./componentApi";
import MsaJob from '../jobs/MsaJob.vue';
import ComponentExceptionsModal from "./ComponentExceptionsModal.vue";

export default defineComponent({
  name: 'NodeContentBody',
  components: {
    ComponentExceptionsModal,
    EsmFoldJob,
    P2RankJob,
    MsaJob,
    draggable,
    DiffDockJob
  },
  props: {
    nodeId: { type: String, required: true },
    name: {
      type: String,
      required: true
    }
  },
  jobsDefinitions: [
    {
      name: "Esmfold light",
      tab: true,
      routeName: "Folding",
      api: componentApi.esmfoldLight
    },
    {
      name: "Esmfold",
      tab: true,
      routeName: "Folding",
      api: componentApi.esmfold
    },
    {
      name: "Rosettafold",
      tab: true,
      routeName: "Folding",
      api: componentApi.rosettafold
    },
    {
      name: "Binding pockets",
      tab: false,
      component: P2RankJob,
      api: componentApi.bindingPockets
    },
    {
      name: "DiffDock",
      tab: false,
      component: DiffDockJob,
      api: componentApi.diffdock
    },
    {
      name: "Msa generation",
      tab: false,
      component: MsaJob,
      api: componentApi.msaGeneration
    },
    {
      name: "Localisation",
      tab: true,
      routeName: "Localisation",
      api: componentApi.localisation
    },
    {
      name: "Solubility",
      tab: true,
      routeName: "Solubility",
      api: componentApi.solubility
    },
    {
      name: "Gene ontology",
      tab: true,
      routeName: "Gene ontology",
      api: componentApi.geneOntology
    },
    {
      name: "Conformations",
      tab: true,
      routeName: "Conformations",
      api: componentApi.conformations
    },
    {
      name: "Protein binder design",
      tab: true,
      routeName: "Protein design",
      api: componentApi.proteinDesign
    },
    {
      name: "Small molecules design",
      tab: true,
      routeName: "Small molecules design",
      api: componentApi.smallMoleculesDesign
    }
  ],
  data() {
    return {
      selectedJobId: '',
      showLastExceptionsModal: false,
      showJobModal: false,
      loading: false,
      error: null as string | null,
      nodeData: null as Node | null,
      jobs: [] as Array<GetJobMetadataResponse & {
        executionStatus: nolabs__application__use_cases__folding__api_models__GetJobStatusResponse | null
      }>,
      results: [] as Array<GetJobMetadataResponse & {
        executionStatus: nolabs__application__use_cases__folding__api_models__GetJobStatusResponse | null
      }>,
      jobErrors: {} as Record<string, string>,
      lastExceptions: [] as string[],
      selectedJobComponent: null as any,
      thumbStyle: {
        right: '4px',
        borderRadius: '7px',
        backgroundColor: '#027be3',
        width: '4px',
        opacity: 0.75
      },
      barStyle: {
        right: '2px',
        borderRadius: '9px',
        backgroundColor: '#027be3',
        width: '8px',
        opacity: 0.2
      }
    };
  },
  async mounted() {
    const workflowStore = useWorkflowStore();
    this.nodeData = workflowStore.getNodeById(this.nodeId);
    await this.updateJobs();
    this.updateErrorsAndExceptions();
    watch(
      () => workflowStore.jobIdsToUpdate,
      async (newJobIds) => {
        for (const jobId of newJobIds) {
          await this.updateJobById(jobId);
        }
      },
      { deep: true }
    );
  },
  watch: {
    'nodeData.data.jobIds': {
      handler: 'updateJobs',
      immediate: true,
      deep: true
    },
    'nodeData.data.last_exceptions': {
      handler: 'updateErrorsAndExceptions',
      immediate: true,
      deep: true
    },
    'nodeData.data.input_property_errors': {
      handle: 'updateErrorsAndExceptions',
      immediate: true,
      deep: true
    },
    'nodeData.data.jobs_errors': {
      handler: 'updateErrorsAndExceptions',
      immediate: true,
      deep: true
    }
  },
  methods: {
    copyContent(s: string) {
      navigator.clipboard.writeText(s);
      Notify.create({
        type: "positive",
        closeBtn: 'Close',
        message: `Copied ${s}`
      });
    },
    async updateJobById(jobId: string) {
      const jobDefinition = this.$options.jobsDefinitions.find(item => item.name === this.name);

      if (!jobDefinition) {
        console.error(`Job definition not found for job: ${this.name}`);
        return;
      }

      // Check if the job is already in the list
      const existingJobIndex = this.jobs.findIndex(job => job.job_id === jobId);
      const existingResultIndex = this.results.findIndex(job => job.job_id === jobId);
      if (existingJobIndex === -1 && existingResultIndex === -1) {
        console.log(`Job ${jobId} is not in the list of jobs for this component.`);
        return;
      }

      try {
        const job = await getJobMetaData(jobId);
        const executionStatus = await jobDefinition.api.executionStatus(jobId);

        // Update the job in this.jobs
        this.jobs[existingJobIndex] = { ...job, executionStatus };

        // Update the results array
        this.jobs = this.jobs.filter(job => job && !job.executionStatus.result_valid) as Array<GetJobMetadataResponse & {
          executionStatus: nolabs__application__use_cases__folding__api_models__GetJobStatusResponse | null
        }>;
        this.results = this.jobs.filter(job => job && job.executionStatus.result_valid) as Array<GetJobMetadataResponse & {
          executionStatus: nolabs__application__use_cases__folding__api_models__GetJobStatusResponse | null
        }>;

        // Update workflow store with running jobs status
        const workflowStore = useWorkflowStore();
        const anyJobRunning = this.jobs.some(job => job.executionStatus?.running);
        if (anyJobRunning) {
          workflowStore.addRunningComponentId(this.nodeId);
        } else {
          workflowStore.removeRunningComponentId(this.nodeId);
        }

        // Remove jobId from the list after updating
        workflowStore.removeJobIdToUpdate(jobId);

      } catch (error) {
        console.error(`Error updating job ${jobId}:`, error);
      }
    },
    async updateJobs() {
      this.loading = true;
      const workflowStore = useWorkflowStore();
      const knownJobIds = new Set(this.jobs.map(job => job.job_id)); // Assuming 'job_id' is the job identifier
      const knownResultIds = new Set(this.results.map(result => result.job_id));

      if (this.nodeData?.data.jobIds) {
        // Fetch new jobs and update existing ones if necessary
        const newJobsWithStatus = await Promise.all(
          this.nodeData.data.jobIds.map(async (jobId: string) => {
            if (!knownJobIds.has(jobId) && !knownResultIds.has(jobId)) {
              const jobDefinition = this.$options.jobsDefinitions.find(item => item.name === this.name);
              const job = await getJobMetaData(jobId);
              const executionStatus = await jobDefinition.api.executionStatus(jobId);
              return { ...job, executionStatus };
            } else {
              return this.jobs.find(job => job.job_id === jobId) || this.results.find(result => result.job_id === jobId);
            }
          })
        );

        // Update the jobs and results lists with the new data
        const newJobs = newJobsWithStatus.filter(job => job && !job.executionStatus.result_valid) as Array<GetJobMetadataResponse & {
          executionStatus: nolabs__application__use_cases__folding__api_models__GetJobStatusResponse | null
        }>;
        const newResults = newJobsWithStatus.filter(job => job && job.executionStatus.result_valid) as Array<GetJobMetadataResponse & {
          executionStatus: nolabs__application__use_cases__folding__api_models__GetJobStatusResponse | null
        }>;

        // Append new jobs and results to existing ones
        this.jobs = [
          ...this.jobs.filter(job => this.nodeData.data.jobIds.includes(job.job_id)),
          ...newJobs.filter(job => !this.jobs.some(existingJob => existingJob.job_id === job.job_id))
        ];

        this.results = [
          ...this.results.filter(result => this.nodeData.data.jobIds.includes(result.job_id)),
          ...newResults.filter(result => !this.results.some(existingResult => existingResult.job_id === result.job_id))
        ];
      }

      this.loading = false;

      // Update workflow store with running jobs status
      const anyJobRunning = this.jobs.some(job => job.executionStatus?.running);
      if (anyJobRunning) {
        workflowStore.addRunningComponentId(this.nodeId);
      } else {
        workflowStore.removeRunningComponentId(this.nodeId);
      }
    },
    async deleteJob(job: GetJobMetadataResponse) {
      this.jobs = this.jobs.filter(j => j.job_id !== job.job_id);
      this.results = this.results.filter(j => j.job_id !== job.job_id);
      const workflowStore = useWorkflowStore();
      await workflowStore.deleteJob(job.job_id);
    },
    updateErrorsAndExceptions() {
      if (this.nodeData?.input_property_errors && this.nodeData.input_property_errors.length > 0) {
        this.lastExceptions = this.nodeData?.input_property_errors.map(x => x.msg + ' - ' + x.loc.join(', ')) || [];
      } else {
        this.lastExceptions = this.nodeData?.data.last_exceptions || [];
      }

      if (this.nodeData?.data.jobs_errors) {
        this.jobErrors = this.nodeData?.data.jobs_errors.reduce((acc: Record<string, string>, error: {
          msg: string;
          job_id: string
        }) => {
          acc[error.job_id] = error.msg;
          return acc;
        }, {}) || {};
      }
    },
    openJob(job: GetJobMetadataResponse) {
      const jobDefinition = this.$options.jobsDefinitions.find(item => item.name === this.name);

      if (jobDefinition) {
        if (jobDefinition.component) {
          this.selectedJobId = job.job_id;
          this.selectedJobComponent = jobDefinition.component;
          this.showJobModal = true;
        } else {
          const routeData = this.$router.resolve({
            name: jobDefinition.routeName, params: { jobId: job.job_id }
          });
          window.open(routeData.href, '_blank');
        }
      }
    },
    updateJobOrder() {
      // Logic to update the job order in the backend if necessary
    },
  },
});
</script>

<style scoped>
.job-section {
  border: 1px dashed #ccc;
  padding: 10px;
}

.result-section {
  border: 1px dashed #ccc;
  padding: 10px;
}

.exception-section {
  border: 1px dashed #ccc;
  padding: 10px;
}

.text-truncate {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}
</style>
