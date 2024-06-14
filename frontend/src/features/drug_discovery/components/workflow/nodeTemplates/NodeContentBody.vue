<template>
  <q-card-section class="job-section q-pa-sm">
    <div class="text-white q-pa-sm">Jobs queue</div>
    <draggable class="q-pa-sm" v-model="jobs" handle=".drag-handle" @end="updateJobOrder" item-key="job_id">
      <template #item="{ element }">
        <q-item class="bg-grey-3 text-black q-pa-sm q-mb-sm q-border-radius-md">
          <q-item-section>
            <q-icon name="drag_indicator" class="drag-handle q-mr-sm" />
          </q-item-section>
          <q-item-section>
            <q-item-label>{{ element.job_id }}</q-item-label>
          </q-item-section>
          <q-item-section>
            <q-item-label>{{ element.job_name }}</q-item-label>
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
            <q-btn @click="openJob(element)" label="View Job" dense />
          </q-item-section>
        </q-item>
      </template>
    </draggable>
  </q-card-section>

  <q-card-section class="result-section q-pa-sm">
    <div class="text-white q-pa-sm">Results</div>
    <draggable class="q-pa-sm" v-model="results" handle=".drag-handle" item-key="job_id">
      <template #item="{ element }">
        <q-item class="bg-grey-3 text-black q-pa-sm q-mb-sm q-border-radius-md">
          <q-item-section>
            <q-icon name="drag_indicator" class="drag-handle q-mr-sm" />
          </q-item-section>
          <q-item-section>
            <q-item-label>{{ element.job_id }}</q-item-label>
          </q-item-section>
          <q-item-section>
            <q-item-label>{{ element.job_name }}</q-item-label>
          </q-item-section>
          <q-item-section>
            <q-item-label>Completed</q-item-label>
          </q-item-section>
          <q-item-section v-if="jobErrors[element.job_id]">
            <q-item-label class="text-red">{{ jobErrors[element.job_id] }}</q-item-label>
          </q-item-section>
          <q-item-section>
            <q-btn @click="openJob(element)" label="View Job" dense />
          </q-item-section>
        </q-item>
      </template>
    </draggable>
  </q-card-section>

  <q-card-section class="exception-section q-pa-sm" v-if="lastExceptions.length">
    <div class="text-white q-pa-sm">Last Exceptions</div>
    <q-item class="bg-grey-3 text-black q-pa-sm q-mb-sm q-border-radius-md" v-for="(exception, index) in lastExceptions"
      :key="index">
      <q-item-section>
        <q-item-label class="text-red">{{ exception }}</q-item-label>
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
import { useWorkflowStore, Node } from 'src/features/drug_discovery/components/workflow/storage';
import { getFoldingJobApi, getFoldingJobStatus, getBindingPocketJobApi, getBindingPocketJobStatus, getDiffDockJobApi, getDiffDockJobStatus } from 'src/features/drug_discovery/refinedApi';
import { GetJobMetadataResponse, nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse } from 'src/refinedApi/client';
import EsmFoldJob from '../jobs/EsmFoldJob.vue';
import P2RankJob from '../jobs/P2RankJob.vue';
import DiffDockJob from '../jobs/DiffDockJob.vue';
import draggable from 'vuedraggable';

export default defineComponent({
  name: 'NodeContentBody',
  components: {
    EsmFoldJob,
    P2RankJob,
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
  computed: {
    jobsFactory() {
      return [
        {
          name: "Esmfold light component",
          tab: false,
          component: EsmFoldJob
        },
        {
          name: "Protein binding pockets prediction",
          tab: false,
          component: P2RankJob
        },
        {
          name: "DiffDock",
          tab: false,
          component: DiffDockJob
        },
        {
          name: "Localisation",
          tab: true,
          routeName: "Localisation"
        },
        {
          name: "Solubility",
          tab: true,
          routeName: "Solubility"
        },
        {
          name: "Gene ontology",
          tab: true,
          routeName: "Gene ontology"
        },
        {
          name: "Conformations",
          tab: true,
          routeName: "Conformations"
        },
        {
          name: "Protein design",
          tab: true,
          routeName: "Protein design"
        },
        {
          name: "Small molecules design",
          tab: true,
          routeName: "Small molecules design"
        }
      ]
    }
  },
  data() {
    return {
      selectedJobId: '',
      showJobModal: false,
      loading: false,
      error: null as string | null,
      nodeData: null as Node | null,
      jobs: [] as Array<GetJobMetadataResponse & { executionStatus: nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse | null }>,
      results: [] as Array<GetJobMetadataResponse & { executionStatus: nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse | null }>,
      jobErrors: {} as Record<string, string>,
      lastExceptions: [] as string[],
      selectedJobComponent: null as any,
    };
  },
  async created() {
    const workflowStore = useWorkflowStore();
    this.nodeData = workflowStore.getNodeById(this.nodeId);
    this.lastExceptions = this.nodeData?.data.last_exceptions || [];
    this.jobErrors = this.nodeData?.data.jobs_errors.reduce((acc: Record<string, string>, error: { msg: string; job_id: string }) => {
      acc[error.job_id] = error.msg;
      return acc;
    }, {}) || {};
    await this.updateJobs();
  },
  watch: {
    'nodeData.data.jobIds': {
      handler: 'updateJobs',
      immediate: true,
      deep: true
    }
  },
  methods: {
    async updateJobs() {
      this.loading = true;
      const diffDockJob = await getDiffDockJobApi("27102be2-6ce7-4fa6-8dd8-c41438cd7012");
      const executionStatus = await getDiffDockJobStatus("27102be2-6ce7-4fa6-8dd8-c41438cd7012");
      let jobsWithStatus;
      if (this.name == 'DiffDock') {
        jobsWithStatus = [{ ...diffDockJob, executionStatus }]
      } else {
        jobsWithStatus = await Promise.all(this.nodeData?.data.jobIds.map(async (jobId: string) => {
          let job;
          let executionStatus;

          switch (this.name) {
            case 'Esmfold light component' || 'Esmfold component' || 'Rosettafold component':
              job = await getFoldingJobApi(jobId);
              executionStatus = await getFoldingJobStatus(jobId);
              break;
            case 'Protein binding pockets prediction':
              job = await getBindingPocketJobApi(jobId);
              executionStatus = await getBindingPocketJobStatus(jobId);
              break;
            // Add cases for other job types and their respective API calls
            // case 'anotherJobType':
            //   job = await getAnotherJobTypeApi(jobId);
            //   executionStatus = await getAnotherJobTypeStatus(jobId);
            //   break;
            default:
              console.error(`Unknown job type: ${this.name}`);
              return null; // Handle unknown job types
          }

          return { ...job, executionStatus };
        }));
      }
      this.jobs = jobsWithStatus.filter(job => job && job.executionStatus && (job.executionStatus.running || !job.result || job.result.length === 0)) as Array<GetJobMetadataResponse & { executionStatus: nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse | null }>;
      this.results = jobsWithStatus.filter(job => job && job.executionStatus && job.result && job.result.length > 0) as Array<GetJobMetadataResponse & { executionStatus: nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse | null }>;

      this.loading = false;
    },
    openJob(job: GetJobMetadataResponse) {
      this.selectedJobId = job.job_id;
      const jobType = this.jobsFactory.find(item => item.name === this.name);

      if (jobType && jobType.component) {
        this.selectedJobComponent = jobType.component;
        this.showJobModal = true;
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
</style>
