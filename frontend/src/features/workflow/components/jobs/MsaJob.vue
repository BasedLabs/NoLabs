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
                  <q-btn v-if="!jobStatus.running" @click="startJob" color="primary" label="Start Job" />
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
          <div class="q-pl-sm q-ma-sm" v-if="jobHasGeneratedData">
            <q-btn @click="downloadA3MFile" color="primary" label="Download .a3m File" />
          </div>
          <div v-else>
            No data available.
          </div>
        </q-card-section>
      </div>
    </div>
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import { QSpinner, QInput, QBtn } from 'quasar';
import {
  nolabs__application__use_cases__msa_generation__api_models__GetJobStatusResponse,
  nolabs__application__use_cases__msa_generation__api_models__JobResponse,
  ProteinContentResponse
} from "../../../../refinedApi/client";
import {changeJobName, getMsaJobApi, getMsajobStatus, getProtein, startMsaJob} from "../../refinedApi";

export default defineComponent({
  name: 'MsaJob',
  props: {
    jobId: String,
  },
  data() {
    return {
      experimentId: null as string | null,
      job: null as nolabs__application__use_cases__msa_generation__api_models__JobResponse | null,
      protein: null as ProteinContentResponse | null,
      jobStatus: null as nolabs__application__use_cases__msa_generation__api_models__GetJobStatusResponse | null,
      editableJobName: '' as string,
    };
  },
  computed: {
    jobHasGeneratedData(): boolean | null {
      return this.job && this.job.result !== null;
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
    this.job = await getMsaJobApi(this.jobId as string);
    if (this.job) {
      this.editableJobName = this.job.job_name || '';
      this.protein = await getProtein(this.job.protein_id);
    }
    this.jobStatus = await getMsajobStatus(this.jobId as string);
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
            message: 'Failed to update job name.'
          });
        }
      }
    },
    async startJob() {
      try {
        await startMsaJob(this.jobId as string);
        // Optionally update job status here if required
      } catch (error) {
        this.$q.notify({
          type: 'negative',
          message: 'Failed to start the job.'
        });
      }
    },
    downloadA3MFile() {
      const blob = new Blob([this.job?.result as string], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${this.protein?.name}.a3m`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    },
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
