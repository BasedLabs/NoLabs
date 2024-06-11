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
          <q-item-section>
            <q-btn @click="openJob(element)" label="View Job" dense />
          </q-item-section>
        </q-item>
      </template>
    </draggable>
  </q-card-section>

  <q-dialog v-model="showJobModal" full-width>
    <q-card style="max-width: 90vw;">
      <esm-fold-job v-if="name === 'Folding'" :job-id="selectedJobId" />
      <!-- Add other job components based on job type -->
      <q-card-actions>
        <q-btn flat label="Close" v-close-popup />
      </q-card-actions>
    </q-card>
  </q-dialog>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import { useWorkflowStore, Node } from 'src/features/drug_discovery/components/workflow/storage';
import { getFoldingJobApi, getFoldingJobStatus } from 'src/features/drug_discovery/refinedApi';
import { GetJobMetadataResponse, nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse } from 'src/refinedApi/client';
import EsmFoldJob from '../jobs/EsmFoldJob.vue';
import draggable from 'vuedraggable';

export default defineComponent({
  name: 'NodeContentBody',
  components: {
    EsmFoldJob,
    draggable,
  },
  props: {
    nodeId: { type: String, required: true },
    name: {
      type: String,
      required: true
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
    };
  },
  async created() {
    const workflowStore = useWorkflowStore();
    this.nodeData = workflowStore.getNodeById(this.nodeId);
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
      const jobsWithStatus = await Promise.all(this.nodeData?.data.jobIds.map(async (jobId: string) => {
        let job;
        let executionStatus;

        switch (this.name) {
          case 'Folding':
            job = await getFoldingJobApi(jobId);
            executionStatus = await getFoldingJobStatus(jobId);
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

      this.jobs = jobsWithStatus.filter(job => job && job.executionStatus && (job.executionStatus.running || !job.result || job.result.length === 0)) as Array<GetJobMetadataResponse & { executionStatus: nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse | null }>;
      this.results = jobsWithStatus.filter(job => job && job.executionStatus && job.result && job.result.length > 0) as Array<GetJobMetadataResponse & { executionStatus: nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse | null }>;

      this.loading = false;
    },
    openJob(job: GetJobMetadataResponse) {
      this.selectedJobId = job.job_id;
      this.showJobModal = true;
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
</style>
