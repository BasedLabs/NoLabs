<template>
  <q-card-section class="job-section q-pa-sm">
    <div class="text-white q-pa-sm text-caption">Jobs queue</div>
    <q-scroll-area v-if="jobs.length > 0" style="height: 30vh">
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

  <q-card-section class="result-section q-pa-sm">
    <div class="text-white q-pa-sm text-caption">Completed jobs</div>
    <q-scroll-area v-if="results.length > 0" style="height: 30vh">
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
    <q-item class="text-black q-pa-sm q-mb-sm q-border-radius-md" v-for="(exception, index) in lastExceptions"
      :key="index">
      <q-item-section>
        <q-btn v-if="exception" color="red" class="q-pm-md">
          <q-icon left name="warning" />
          <div>{{ exception }}</div>
        </q-btn>
      </q-item-section>
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
import { defineComponent } from 'vue';
import { useWorkflowStore, Node } from '../storage';
import {
  GetJobMetadataResponse,
  nolabs__application__use_cases__folding__api_models__GetJobStatusResponse
} from 'src/refinedApi/client';
import EsmFoldJob from '../jobs/EsmFoldJob.vue';
import P2RankJob from '../jobs/P2RankJob.vue';
import DiffDockJob from '../jobs/DiffDockJob.vue';
import draggable from 'vuedraggable';
import { Notify } from "quasar";
import componentApi from "./componentApi";
import MsaJob from '../jobs/MsaJob.vue';

export default defineComponent({
  name: 'NodeContentBody',
  components: {
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
      name: "Protein design",
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
    };
  },
  async mounted() {
    const workflowStore = useWorkflowStore();
    this.nodeData = workflowStore.getNodeById(this.nodeId);
    await this.updateJobs();
    this.updateErrorsAndExceptions();
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
    async updateJobs() {
      this.loading = true;
      const workflowStore = useWorkflowStore();
      let jobsWithStatus = [];
      if (this.nodeData?.data.jobIds) {
        jobsWithStatus = await Promise.all(this.nodeData?.data?.jobIds?.map(async (jobId: string) => {
          let job;
          let executionStatus;

          const jobDefinition = this.$options.jobsDefinitions.find(item => item.name === this.name);

          job = await jobDefinition.api.getJob(jobId);
          executionStatus = await jobDefinition.api.executionStatus(jobId);

          return { ...job, executionStatus };
        }));
      }
      this.jobs = jobsWithStatus.filter(job => job && !job.executionStatus.result_valid) as Array<GetJobMetadataResponse & {
        executionStatus: nolabs__application__use_cases__folding__api_models__GetJobStatusResponse | null
      }>;
      this.results = jobsWithStatus.filter(job => job && job.executionStatus.result_valid) as Array<GetJobMetadataResponse & {
        executionStatus: nolabs__application__use_cases__folding__api_models__GetJobStatusResponse | null
      }>;

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
      const workflowStore = useWorkflowStore();
      await workflowStore.deleteJob(job.job_id);
      this.updateJobs();
    },
    updateErrorsAndExceptions() {
      this.lastExceptions = this.nodeData?.data.last_exceptions || [];
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
          this.$router.push({
            name: jobDefinition.routeName, params: { jobId: job.job_id }
          });
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
