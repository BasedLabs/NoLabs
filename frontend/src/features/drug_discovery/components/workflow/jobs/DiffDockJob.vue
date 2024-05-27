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
              <q-item-section>{{ job?.job_name }}</q-item-section>
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
            <PdbViewer :pdb-file="protein?.pdb_content" :key="protein?.name"/>
          </div>
        </q-card-section>
      </div>
    </div>
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import { QSpinnerOrbit, QSpinner } from 'quasar';
import PdbViewer from 'src/components/PdbViewer.vue';
import {
  nolabs__refined__application__use_cases__folding__api_models__JobResponse, ProteinResponse,
  nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse
} from "../../../../../refinedApi/client";
import { getFoldingJobApi, getProtein, getFoldingJobStatus } from "../../../refinedApi";

export default defineComponent({
  name: 'DiffDockJob',
  props: {
    jobId: String,
  },
  data() {
    return {
      job: null as nolabs__refined__application__use_cases__folding__api_models__JobResponse | null,
      protein: null as ProteinResponse | null,
      jobStatus: null as nolabs__refined__application__use_cases__folding__api_models__GetJobStatusResponse | null,
    };
  },
  computed: {
    jobHasGeneratedData(): boolean {
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
    this.$q.loading.show({
      spinner: QSpinnerOrbit,
      message: `Loading Experiment ${this.jobId}`,
    });

    this.job = await getFoldingJobApi(this.jobId as string);
    if (this.job && this.job.proteins.length > 0) {
      this.protein = await getProtein(this.job.proteins[0]);
    }

    this.jobStatus = await getFoldingJobStatus(this.jobId as string);

    this.$q.loading.hide();
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
