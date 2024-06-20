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
            <q-item>
              <q-item-section>
                <q-item-label>Folding Backend</q-item-label>
              </q-item-section>
              <q-item-section>{{ job?.backend }}</q-item-section>
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
            <PdbViewer :pdb-file="pdbFile" :key="protein?.name"/>
          </div>
        </q-card-section>
      </div>
    </div>
  </q-card>
</template>


<script lang="ts">
import { defineComponent } from 'vue';
import { QSpinnerOrbit, QSpinner, QInput } from 'quasar';
import PdbViewer from 'src/components/PdbViewer.vue';
import {
  nolabs__application__use_cases__folding__api_models__JobResponse, ProteinContentResponse,
  nolabs__application__use_cases__folding__api_models__GetJobStatusResponse,
  nolabs__application__use_cases__folding__api_models__SetupJobRequest
} from "src/refinedApi/client";
import { getFoldingJobApi, getProteinContent, getFoldingJobStatus, setupFoldingJob, changeJobName } from "src/features/workflow/refinedApi";

export default defineComponent({
  name: 'EsmFoldJob',
  props: {
    jobId: String,
  },
  data() {
    return {
      experimentId: null as string | null,
      job: null as nolabs__application__use_cases__folding__api_models__JobResponse | null,
      protein: null as ProteinContentResponse | null,
      jobStatus: null as nolabs__application__use_cases__folding__api_models__GetJobStatusResponse | null,
      editableJobName: '' as string,
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
    pdbFile(): File {
      if (this.jobHasGeneratedData && this.job?.result[0].pdb) {
        return new File([new Blob([this.job?.result[0].pdb])], this.protein?.name + ".pdb");
      }
      return new File([], "empty.pdb");
    }
  },
  async mounted() {

    this.experimentId = this.$route.params.experimentId as string;

    this.job = await getFoldingJobApi(this.jobId as string);
    if (this.job) {
      this.editableJobName = this.job.job_name || '';
    }
    if (this.job && this.job.protein_ids.length > 0) {
      this.protein = await getProteinContent(this.job.protein_ids[0]);
    }

    this.jobStatus = await getFoldingJobStatus(this.jobId as string);
  },
  methods: {
    async updateJobName() {
      if (this.job) {
        try {
          await changeJobName(this.job.job_id, this.editableJobName);        }
        catch (error) {
          this.editableJobName = this.job.job_name;
          this.$q.notify({
            type: 'negative',
            message: 'Failed to update job name.'
          });
        }
      }
    },
    async updateJobParams() {
      if (this.job) {
        const setupJobRequest: nolabs__application__use_cases__folding__api_models__SetupJobRequest = {
          experiment_id: this.experimentId as string,
          backend: this.job.backend,
          protein_ids: this.job.protein_ids,
          job_id: this.job.job_id,
          job_name: this.editableJobName,
        };

        try {
          this.job = await setupFoldingJob(setupJobRequest);
          this.$q.notify({
            type: 'positive',
            message: 'Job name updated successfully!'
          });
        } catch (error) {
          this.$q.notify({
            type: 'negative',
            message: 'Failed to update job name.'
          });
        }
      }
    }
  },
  components: {
    PdbViewer
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
