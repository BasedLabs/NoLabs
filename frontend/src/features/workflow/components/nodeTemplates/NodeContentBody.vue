<template>
  <q-card-section class="job-section q-pa-sm" style="width: 100%;">
    <div class="text-white q-pa-sm text-caption">Jobs queue</div>
    <q-scroll-area
      v-if="jobs.length > 0"
      visible
      :thumb-style="thumbStyle"
      :bar-style="barStyle"
      style="height: 30vh; width: 100%;"
    >
      <draggable
        class="q-pa-sm"
        v-model="jobs"
        handle=".drag-handle"
        @end="updateJobOrder"
        item-key="job_id"
      >
        <template #item="{ element }">
          <q-item class="bg-grey-3 text-black">
            <!-- Copy ID Button -->
            <q-item-section style="flex: 0 0 60px;">
              <q-btn
                @click="copyContent(element.job_id)"
                icon="content_copy"
                label="ID"
                dense
              />
            </q-item-section>

            <!-- Job Name -->
            <q-item-section class="col">
              <q-item-label class="job-name">{{ element.job_name }}</q-item-label>
              <q-tooltip>{{ element.job_name }}</q-tooltip>
            </q-item-section>

            <!-- Execution Status -->
            <q-item-section style="flex: 0 0 120px;">
              <div v-if="element.executionStatus === null">No status</div>
              <div v-else>
                <!-- Display Error for FAILED state -->
                <div v-if="element.executionStatus.state === 'FAILED'">
                  <q-btn dense label="ERROR" text-color="black" icon="warning" color="yellow"></q-btn>
                  <q-tooltip>{{ element.executionStatus.state_message }}</q-tooltip>
                </div>
                <!-- Display Running... for RUNNING state -->
                <div v-else-if="element.executionStatus.state === 'RUNNING'">
                  <q-spinner color="primary" size="20px" />
                  Running...
                </div>
                <!-- Handle other states -->
                <div v-else>
                  {{ element.executionStatus.state }}
                </div>
              </div>
            </q-item-section>

            <!-- Job Error Message -->
            <q-item-section
              v-if="jobErrors[element.job_id]"
              style="flex: 0 0 150px;"
            >
              <q-item-label class="text-red">{{ jobErrors[element.job_id] }}</q-item-label>
            </q-item-section>

            <!-- View Button -->
            <q-item-section style="flex: 0 0 60px;">
              <q-btn @click="openJob(element)" label="View" dense />
            </q-item-section>

            <!-- Delete Button -->
            <q-item-section style="flex: 0 0 70px;">
              <q-btn
                @click="deleteJob(element)"
                color="negative"
                label="Delete"
                dense
              />
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

  <!-- Component Error Section -->
  <q-card-section class="component-error-section q-pa-sm" v-if="nodeData?.data.state === 'FAILED'">
    <div class="text-white q-pa-sm">Component Error</div>
    <q-item class="text-black q-pa-sm q-mb-sm q-border-radius-md">
      <q-item-section>
        <q-btn color="yellow" text-color="black" label="Component Error" class="q-pm-md">
          <q-icon left name="warning" />
        </q-btn>
        <div class="error-message">
          {{ nodeData.data.stateMessage }}
        </div>
      </q-item-section>
    </q-item>
  </q-card-section>

  <!-- Job Modal -->
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
import { useWorkflowStore, Node } from 'src/features/workflow/components/storage';
import {
  getJobMetaData,
} from 'src/features/workflow/refinedApi';
import {
  GetJobMetadataResponse,
  JobStateEnum,
  ComponentStateEnum, // Import ComponentStateEnum
} from 'src/refinedApi/client';
import EsmFoldJob from '../jobs/EsmFoldJob.vue';
import P2RankJob from '../jobs/P2RankJob.vue';
import DiffDockJob from '../jobs/DiffDockJob.vue';
import draggable from 'vuedraggable';
import { Notify } from 'quasar';
import componentApi from './componentApi';
import MsaJob from '../jobs/MsaJob.vue';
import ComponentExceptionsModal from './ComponentExceptionsModal.vue';

export default defineComponent({
  name: 'NodeContentBody',
  components: {
    ComponentExceptionsModal,
    EsmFoldJob,
    P2RankJob,
    MsaJob,
    draggable,
    DiffDockJob,
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
      name: "Blast",
      tab: true,
      routeName: "Blast",
      api: componentApi.blast
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
    },
    {
      name: "ProteinMPNN",
      tab: true,
      routeName: "ProteinMPNN",
      api: componentApi.proteinMpnn
    }
  ],
  data() {
    return {
      error: null as string | null,
      nodeData: null as Node | null,
      jobs: [] as Array<GetJobMetadataResponse & {
        executionStatus: JobStateEnum | null
      }>,
      results: [] as Array<GetJobMetadataResponse & {
        executionStatus: JobStateEnum | null
      }>,
      jobErrors: {} as Record<string, string>,
      selectedJobComponent: null as any,
      thumbStyle: {
        right: '4px',
        borderRadius: '7px',
        backgroundColor: '#027be3',
        width: '4px',
        opacity: 0.75,
      },
      barStyle: {
        right: '2px',
        borderRadius: '9px',
        backgroundColor: '#027be3',
        width: '8px',
        opacity: 0.2,
      },
      JobStateEnum, // Make JobStateEnum available in the template
      showJobModal: false,
      selectedJobId: '',
      ComponentStateEnum, // Make ComponentStateEnum available in the template
    };
  },
  async mounted() {
    const workflowStore = useWorkflowStore();
    this.nodeData = workflowStore.getNodeById(this.nodeId);
    await this.updateJobs();
  },
  watch: {
    'nodeData.data.jobIdsToUpdate': {
      handler: 'updateJobsToUpdate',
      immediate: true,
      deep: true
    },
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
    async updateJobsToUpdate() {
      if (!this.nodeData || !this.nodeData.data.jobIdsToUpdate || this.nodeData.data.jobIdsToUpdate.length === 0) {
        return;
      }

      // Loop over each job ID in the jobIdsToUpdate and call updateJobById for each
      for (const jobId of this.nodeData.data.jobIdsToUpdate) {
        await this.updateJobById(jobId);
      }
    },
    async updateJobById(jobId: string) {
      const jobDefinition = this.$options.jobsDefinitions.find(item => item.name === this.name);

      if (!jobDefinition) {
        console.error(`Job definition not found for job: ${this.name}`);
        return;
      }

      try {
        const job = await getJobMetaData(jobId);
        const executionStatus = await jobDefinition.api.executionStatus(jobId);

        // Check if the job is already in the jobs or results lists
        const existingJobIndex = this.jobs.findIndex(job => job.job_id === jobId);
        const existingResultIndex = this.results.findIndex(job => job.job_id === jobId);

        // If the job is completed, it should be in the results list and removed from jobs list
        if (executionStatus.state === JobStateEnum.COMPLETED) {
          if (existingJobIndex !== -1) {
            // Remove the job from the jobs list if it's there
            this.jobs.splice(existingJobIndex, 1);
          }
          if (existingResultIndex === -1) {
            // If it's not already in results, add it
            this.results.push({ ...job, executionStatus });
          } else {
            // Update the job in results if it's already there
            this.results[existingResultIndex] = { ...job, executionStatus };
          }
        } else {
          // If the job is not completed, it should be in the jobs list and removed from results list
          if (existingResultIndex !== -1) {
            // Remove the job from the results list if it's there
            this.results.splice(existingResultIndex, 1);
          }
          if (existingJobIndex === -1) {
            // If it's not already in jobs, add it
            this.jobs.push({ ...job, executionStatus });
          } else {
            // Update the job in jobs if it's already there
            this.jobs[existingJobIndex] = { ...job, executionStatus };
          }
        }

        // Update workflow store with running jobs status
        const workflowStore = useWorkflowStore();
        // Remove jobId from the list after updating
        workflowStore.removeJobIdToUpdate(this.nodeId, jobId);

      } catch (error) {
        console.error(`Error updating job ${jobId}:`, error);
      }
    },
    async updateJobs() {
      this.loading = true;

      if (this.nodeData?.data.jobIds) {
        // Fetch new jobs and update existing ones if necessary
        const newJobsWithStatus = await Promise.all(
          this.nodeData.data.jobIds.map(async (jobId: string) => {
              const jobDefinition = this.$options.jobsDefinitions.find(item => item.name === this.name);
              const job = await getJobMetaData(jobId);
              const executionStatus = await jobDefinition.api.executionStatus(jobId);
              return { ...job, executionStatus };
          })
        );

        // Update the jobs and results lists with the new data
        const newJobs = newJobsWithStatus.filter(job => job && job.executionStatus?.state !== JobStateEnum.COMPLETED) as Array<GetJobMetadataResponse & {
          executionStatus: JobStateEnum | null;
        }>;
        const newResults = newJobsWithStatus.filter(job => job && job.executionStatus?.state === JobStateEnum.COMPLETED) as Array<GetJobMetadataResponse & {
          executionStatus: JobStateEnum | null;
        }>;

        // Append new jobs and results to existing ones
        this.jobs = [
          ...this.jobs.filter(job => this.nodeData.data.jobIds.includes(job.job_id)),
          ...newJobs.filter(job => !this.jobs.some(existingJob => existingJob.job_id === job.job_id)),
        ];

        this.results = [
          ...this.results.filter(result => this.nodeData.data.jobIds.includes(result.job_id)),
          ...newResults.filter(result => !this.results.some(existingResult => existingResult.job_id === result.job_id)),
        ];
      }

      this.loading = false;
    },
    async deleteJob(job: GetJobMetadataResponse) {
      this.jobs = this.jobs.filter(j => j.job_id !== job.job_id);
      this.results = this.results.filter(j => j.job_id !== job.job_id);
      const workflowStore = useWorkflowStore();
      await workflowStore.deleteJob(job.job_id, this.nodeId);
    },
    openJob(job: GetJobMetadataResponse) {
      const jobDefinition = this.$options.jobsDefinitions.find(item => item.name === this.name);

      if (jobDefinition) {
        if (jobDefinition.component) {
          this.selectedJobId = job.job_id;
          this.selectedJobComponent = jobDefinition.component;
          this.showJobModal = true;
        } else {
          const experimentId = this.$route.params.experimentId as string;
          const routeData = this.$router.resolve({
            name: jobDefinition.routeName,
            params: { experimentId: experimentId, jobId: job.job_id },
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

.component-error-section {
  border: 1px dashed #ccc;
  padding: 10px;
}

.job-name {
  display: -webkit-box;
  -webkit-line-clamp: 2; /* Limit to 2 lines */
  -webkit-box-orient: vertical;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: normal; /* Allow text to wrap */
}

.error-message {
  margin-top: 10px;
  color: red;
}
</style>
