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
                <q-input v-model="editableJobName" @keyup.enter="updateJobName" dense clearable />
                <q-btn v-if="editableJobName !== (job?.job_name || '')" icon="check" color="info" flat
                  @click="updateJobName" />
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
                  <q-btn v-if="!jobStatus.running" @click="startJob" color="info" label="Start" />
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
                <q-item-label>Descriptions</q-item-label>
              </q-item-section>
              <q-item-section>{{ job?.descriptions }}</q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Alignments</q-item-label>
              </q-item-section>
              <q-item-section>{{ job?.alignments }}</q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Hitlist Size</q-item-label>
              </q-item-section>
              <q-item-section>{{ job?.hitlist_size }}</q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Expect</q-item-label>
              </q-item-section>
              <q-item-section>{{ job?.expect }}</q-item-section>
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
          <div v-if="jobHasGeneratedData">
            <div v-for="result in job?.result" :key="result.protein_id">
              <p>Program: {{ result.program }}</p>
              <p>Database: {{ result.database }}</p>
              <p>Query ID: {{ result.query_id }}</p>
              <p>Query Definition: {{ result.query_def }}</p>
              <p>Query Length: {{ result.query_len }}</p>
              <div v-for="hit in result.hits" :key="hit.id" class="hit-section">
                <h6>Hit ID: {{ hit.id }}</h6>
                <p>Definition: {{ hit.definition }}</p>
                <p>Accession: <a :href="'https://www.rcsb.org/structure/' + hit.accession" target="_blank">{{ hit.accession }}</a></p>
                <p>Length: {{ hit.length }}</p>
                <div v-for="hsp in hit.hsps" :key="hsp.num" class="hsp-section">
                  <p>HSP Number: {{ hsp.num }}</p>
                  <p>Bit Score: {{ hsp.bit_score }}</p>
                  <p>Score: {{ hsp.score }}</p>
                  <p>E-value: {{ hsp.evalue }}</p>
                  <p>Query From: {{ hsp.query_from }}</p>
                  <p>Query To: {{ hsp.query_to }}</p>
                  <p>Hit From: {{ hsp.hit_from }}</p>
                  <p>Hit To: {{ hsp.hit_to }}</p>
                  <p>Query Frame: {{ hsp.query_frame }}</p>
                  <p>Hit Frame: {{ hsp.hit_frame }}</p>
                  <p>Identity: {{ hsp.identity }}</p>
                  <p>Positive: {{ hsp.positive }}</p>
                  <p>Gaps: {{ hsp.gaps }}</p>
                  <p>Alignment Length: {{ hsp.align_len }}</p>
                  <p>Query Sequence: <span class="sequence">{{ hsp.qseq }}</span></p>
                  <p>Hit Sequence: <span class="sequence">{{ hsp.hseq }}</span></p>
                  <p>Midline: <span class="sequence">{{ hsp.midline }}</span></p>
                </div>
              </div>
            </div>
          </div>
          <div v-else>
            No results available.
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
  nolabs__application__use_cases__blast__api_models__JobResponse, 
  ProteinContentResponse, 
  nolabs__application__use_cases__blast__api_models__GetJobStatusResponse, 
  nolabs__application__use_cases__blast__api_models__SetupJobRequest 
} from "src/refinedApi/client";
import { 
  getBlastJobApi, 
  getProteinContent, 
  getBlastJobStatus, 
  changeJobName, 
  setupBlastJob, 
  startBlastJob 
} from "src/features/workflow/refinedApi";

export default defineComponent({
  name: 'BlastJob',
  props: {
    jobId: String,
  },
  data() {
    return {
      experimentId: null as string | null,
      job: null as nolabs__application__use_cases__blast__api_models__JobResponse | null,
      protein: null as ProteinContentResponse | null,
      jobStatus: null as nolabs__application__use_cases__blast__api_models__GetJobStatusResponse | null,
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
  },
  async mounted() {
    this.experimentId = this.$route.params.experimentId as string;

    this.job = await getBlastJobApi(this.jobId as string);
    if (this.job) {
      this.editableJobName = this.job.job_name || '';
    }

    this.protein = await getProteinContent(this.job?.protein_id);

    this.jobStatus = await getBlastJobStatus(this.jobId as string);
  },
  methods: {
    async updateJobName() {
      if (this.job && this.editableJobName !== this.job.job_name) {
        try {
          await changeJobName(this.job.job_id, this.editableJobName);
          this.$q.notify({
            type: 'positive',
            message: 'Job name updated successfully.',
          });
        } catch (error) {
          this.editableJobName = this.job.job_name;
          this.$q.notify({
            type: 'negative',
            message: 'Failed to update job name.',
          });
        }
      }
    },
    async startJob() {
      try {
        await startBlastJob(this.jobId as string);
        this.$q.notify({
          type: 'positive',
          message: 'Job started successfully.',
        });
        this.jobStatus = await getBlastJobStatus(this.jobId as string);
      } catch (error) {
        this.$q.notify({
          type: 'negative',
          message: 'Failed to start the job.',
        });
      }
    }
  },
  components: {
    QBtn, // Ensure QBtn component is registered
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

.hit-section {
  border: 1px solid #00e1ff;
  padding: 10px;
  margin-bottom: 10px;
}

.hsp-section {
  border: 1px solid #00e1ff;
  padding: 10px;
  margin-bottom: 10px;
}

.sequence {
  font-family: monospace;
  background-color: #000000;
  padding: 2px 4px;
  border-radius: 4px;
}
</style>
