<template>
    <q-card-section class="job-section q-pa-sm">
      <h7 class="text-white q-pa-sm">Jobs queue</h7> 
      <q-btn class="q-pa-sm" @click="addJob" label="+ Add Job" color="black" />
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
            </q-item-section>
            <q-item-section>
              <q-item-label>{{ element.job_name }}</q-item-label>
            </q-item-section>
            <q-item-section>
              <div v-if="element.executionStatus === null">Loading...</div>
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
  
    <q-card-section>
      <h7 class="text-white">Results</h7>
      <q-list v-if="results && results.length">
        <q-item v-for="result in results" :key="result.id">
          <q-item-section>{{ result.name }}</q-item-section>
        </q-item>
      </q-list>
    </q-card-section>
  
    <q-dialog v-model="showJobModal" full-width>
      <q-card style="max-width: 90vw;">
        <esm-fold-job :job-id="selectedJobId" />
        <q-card-actions align="right">
          <q-btn flat label="Close" v-close-popup />
        </q-card-actions>
      </q-card>
    </q-dialog>
  </template>
  
  <script lang="ts">
  import { defineComponent } from 'vue';
  import { getFoldingJobsApi, getFoldingJobStatus } from "../../../refinedApi";
  import { GetJobMetadataResponse, nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse } from "../../../../refinedApi/client";
  import EsmFoldJob from "../jobs/EsmFoldJob.vue";
  import draggable from 'vuedraggable';
  
  export default defineComponent({
    name: 'NodeContentBody',
    components: {
      EsmFoldJob,
      draggable,
    },
    props: {
      results: Array,
    },
    data() {
      return {
        jobs: [] as Array<GetJobMetadataResponse & { executionStatus: nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse | null }>,
        selectedJobId: '',
        showJobModal: false,
        loading: false,
      };
    },
    async mounted() {
      this.loading = true;
      const experimentId = this.$route.params.experimentId as string;
      const jobs = await getFoldingJobsApi(experimentId);
  
      // Fetch execution status for each job
      const jobsWithStatus = await Promise.all(jobs.map(async job => {
        const executionStatus = await getFoldingJobStatus(job.job_id);
        return { ...job, executionStatus };
      }));
  
      this.jobs = jobsWithStatus;
      this.loading = false;
    },
    methods: {
      addJob() {
        // Logic to add a new job
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
  </style>
  