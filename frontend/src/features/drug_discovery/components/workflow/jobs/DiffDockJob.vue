<template>
  <q-separator></q-separator>
  <q-card>
    <div class="row q-gutter-md">
      <div class="col-12">
        <q-separator></q-separator>
        <q-card-section>
          <div class="text-h6">Job Details</div>
        </q-card-section>
        <q-card-section>
          <q-list>
            <q-item>
              <q-item-section>
                <q-item-label>Job ID</q-item-label>
              </q-item-section>
              <q-item-section>{{ jobId }}</q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Job Name</q-item-label>
              </q-item-section>
              <q-item-section>
                <q-input v-model="editableJobName" @blur="updateJobName" />
              </q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Job Status</q-item-label>
              </q-item-section>
              <q-item-section>
                <template v-if="jobStatus === null">
                  Loading...
                </template>
                <template v-else>
                  <q-spinner v-if="jobStatus.running" color="primary" size="20px" />
                  {{ jobStatusText }}
                </template>
              </q-item-section>
            </q-item>
          </q-list>
        </q-card-section>
      </div>
      <q-separator></q-separator>
      <div class="col-12">
        <q-separator></q-separator>
        <q-card-section>
          <div class="text-h6">Job Inputs</div>
        </q-card-section>
        <q-card-section>
          <q-list>
            <q-item>
              <q-item-section>
                <q-item-label>Protein FASTA Content</q-item-label>
              </q-item-section>
              <q-item-section class="fasta-content-container">
                <div class="fasta-content">{{ protein?.fasta_content }}</div>
              </q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Number of samples</q-item-label>
              </q-item-section>
              <q-item-section>
                <q-input v-model.number="samplesPerComplex" type="number" @blur="updateNumberOfSamples" />
              </q-item-section>
            </q-item>
          </q-list>
        </q-card-section>
      </div>
      <q-separator></q-separator>
      <div class="col-12">
        <q-separator></q-separator>
        <q-card-section>
          <div class="text-h6">Result</div>
        </q-card-section>
        <q-card-section>
          <DiffDockResult v-if="jobHasGeneratedData" :job="job" :protein="protein" />
        </q-card-section>
      </div>
    </div>
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import { QSpinner, QInput } from 'quasar';
import DiffDockResult from './DiffDockResult.vue';
import {
  nolabs__refined__application__use_cases__diffdock__api_models__JobResponse,
  ProteinResponse,
  nolabs__refined__application__use_cases__diffdock__api_models__GetJobStatusResponse,
  nolabs__refined__application__use_cases__diffdock__api_models__SetupJobRequest
} from "../../../../../refinedApi/client";
import { getDiffDockJobApi, getProtein, getDiffDockJobStatus, changeJobName, setupDiffDockJob } from "../../../refinedApi";

export default defineComponent({
  name: 'DiffDockJob',
  props: {
    jobId: String,
  },
  data() {
    return {
      experimentId: null as string | null,
      job: null as nolabs__refined__application__use_cases__diffdock__api_models__JobResponse | null,
      protein: null as ProteinResponse | null,
      jobStatus: null as nolabs__refined__application__use_cases__diffdock__api_models__GetJobStatusResponse | null,
      editableJobName: '' as string,
      samplesPerComplex: null as number | null
    };
  },
  computed: {
    jobHasGeneratedData(): boolean | null {
      return this.job && this.job.result.length > 0;
    },
    jobStatusText(): string {
      if (this.jobStatus === null) {
        return '';
      }
      return this.jobStatus.running ? 'Running...' : 'Not running';
    },
  },
  async mounted() {
    this.experimentId = this.$route.params.experimentId as string;

    this.job = await getDiffDockJobApi(this.jobId as string);
    if (this.job) {
      this.editableJobName = this.job.job_name || '';
      this.samplesPerComplex = this.job.samples_per_complex || null;
    }

    this.protein = await getProtein(this.job.protein_id);

    this.jobStatus = await getDiffDockJobStatus(this.jobId as string);
  },
  methods: {
    async updateJobName() {
      if (this.job) {
        try {
          await changeJobName(this.job.job_id, this.editableJobName);
        } catch (error) {
          this.editableJobName = this.job.job_name;
          this.$q.notify({
            type: 'negative',
            message: 'Failed to update job name.',
          });
        }
      }
    },
    async updateNumberOfSamples() {
      if (this.job && this.samplesPerComplex !== null) {
        const request: nolabs__refined__application__use_cases__diffdock__api_models__SetupJobRequest = {
          experiment_id: this.experimentId as string,
          protein_id: this.job.protein_id,
          ligand_id: this.job.ligand_id,
          samples_per_complex: this.samplesPerComplex,
          job_id: this.job.job_id,
          job_name: this.job.job_name,
        };
        try {
          await setupDiffDockJob(request);
        } catch (error) {
          this.$q.notify({
            type: 'negative',
            message: 'Failed to update number of samples.',
          });
        }
      }
    }
  },
  components: {
    DiffDockResult,
  },
});
</script>

<style scoped>
.fasta-content-container {
  width: 100%;
  overflow: hidden;
}

.fasta-content {
  white-space: pre-wrap;
  word-wrap: break-word;
  word-break: break-all;
}
</style>
