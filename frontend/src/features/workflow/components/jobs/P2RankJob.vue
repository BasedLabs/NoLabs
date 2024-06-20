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
            <PdbViewer :pdb-file="pdbFile" :pocket-ids="job?.result" :key="protein?.name" />
          </div>
        </q-card-section>
      </div>
    </div>
  </q-card>
</template>


<script lang="ts">
import { defineComponent } from 'vue';
import { QSpinner, QInput } from 'quasar';
import PdbViewer from 'src/components/PdbViewer.vue';
import {
  nolabs__application__use_cases__binding_pockets__api_models__GetJobStatusResponse,
  nolabs__application__use_cases__binding_pockets__api_models__JobResponse,
  ProteinContentResponse
} from "../../../../refinedApi/client";
import {changeJobName, getBindingPocketJobApi, getBindingPocketJobStatus, getProtein} from "../../refinedApi";

export default defineComponent({
  name: 'P2RankJob',
  props: {
    jobId: String,
  },
  data() {
    return {
      experimentId: null as string | null,
      job: null as nolabs__application__use_cases__binding_pockets__api_models__JobResponse | null,
      protein: null as ProteinContentResponse | null,
      jobStatus: null as nolabs__application__use_cases__binding_pockets__api_models__GetJobStatusResponse | null,
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
      if (this.jobHasGeneratedData && this.protein?.pdb_content) {
        return new File([new Blob([this.protein?.pdb_content])], this.protein?.name + ".pdb");
      }
      return new File([], "empty.pdb");
    }
  },
  async mounted() {

    this.experimentId = this.$route.params.experimentId as string;

    this.job = await getBindingPocketJobApi(this.jobId as string);
    if (this.job) {
      this.editableJobName = this.job.job_name || '';
    }

    this.protein = await getProtein(this.job.protein_id);

    this.jobStatus = await getBindingPocketJobStatus(this.jobId as string);

  },
  methods: {
    async updateJobName() {
      if (this.job) {
        try {
          await changeJobName(this.job.job_id, this.editableJobName);
        }
        catch (error) {
          this.editableJobName = this.job.job_name;
          this.$q.notify({
            type: 'negative',
            message: 'Failed to update job name.'
          });
        }
      }
    },
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
