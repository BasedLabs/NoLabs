<template>
  <q-separator></q-separator>
  <q-card>
    <div class="row q-gutter-md">
      <!-- Job Details Section -->
      <div class="col-12">
        <q-separator></q-separator>
        <q-card-section>
          <div class="text-h6">Job Details</div>
        </q-card-section>
        <q-card-section>
          <q-list>
            <!-- Job ID -->
            <q-item>
              <q-item-section>
                <q-item-label>Job ID</q-item-label>
              </q-item-section>
              <q-item-section>{{ jobId }}</q-item-section>
            </q-item>
            <!-- Job Name -->
            <q-item>
              <q-item-section>
                <q-item-label>Job Name</q-item-label>
              </q-item-section>
              <q-item-section>
                <q-input v-model="editableJobName" @keyup.enter="updateJobName" dense clearable />
                <q-btn v-if="editableJobName !== (job?.job_name || '')" icon="check" color="info" flat @click="updateJobName" />
              </q-item-section>
            </q-item>
            <!-- Job Status -->
            <q-item>
              <q-item-section>
                <q-item-label>Job Status</q-item-label>
              </q-item-section>
              <q-item-section>
                <template v-if="jobStatus === null">
                  Loading...
                </template>
                <template v-else>
                  <q-spinner v-if="jobStatus.state === 'RUNNING'" color="primary" size="20px" />
                  {{ jobStatus.state }}
                  <q-btn v-if="jobStatus.state !== 'RUNNING'" @click="startJob" color="info" label="Start" />
                </template>
              </q-item-section>
            </q-item>
          </q-list>
        </q-card-section>
      </div>
      <q-separator></q-separator>
      <!-- Job Inputs Section -->
      <div class="col-12">
        <q-separator></q-separator>
        <q-card-section>
          <div class="text-h6">Job Inputs</div>
        </q-card-section>
        <q-card-section>
          <q-list>
            <!-- Temperature -->
            <q-item>
              <q-item-section>
                <q-item-label>Sampling Temperature</q-item-label>
              </q-item-section>
              <q-item-section>
                <q-input v-model.number="samplingTemp" type="number" @keyup.enter="updateSamplingTemp" dense clearable />
                <q-btn v-if="samplingTemp !== (job?.sampling_temp || null)" icon="check" color="info" flat @click="updateSamplingTemp" />
              </q-item-section>
            </q-item>
            <!-- Number of Sequences per Target -->
            <q-item>
              <q-item-section>
                <q-item-label>Number of Sequences per Target</q-item-label>
              </q-item-section>
              <q-item-section>
                <q-input v-model.number="numSeqPerTarget" type="number" @keyup.enter="updateNumSeqPerTarget" dense clearable />
                <q-btn v-if="numSeqPerTarget !== (job?.num_seq_per_target || null)" icon="check" color="info" flat @click="updateNumSeqPerTarget" />
              </q-item-section>
            </q-item>
            <!-- Is Homomer -->
            <q-item>
              <q-item-section>
                <q-item-label>Is Homomer</q-item-label>
              </q-item-section>
              <q-item-section>
                <q-toggle v-model="isHomomer" @change="updateIsHomomer" />
              </q-item-section>
            </q-item>
          </q-list>
        </q-card-section>
      </div>
      <q-separator></q-separator>
      <!-- Result Section -->
      <div class="col-12">
        <q-separator></q-separator>
        <q-card-section>
          <div class="text-h6">Result</div>
        </q-card-section>
        <q-card-section>
          <ProteinMPNNResult v-if="jobHasGeneratedData" :job="job" :protein="protein" />
        </q-card-section>
      </div>
    </div>
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import { QBtn, QToggle } from 'quasar';
import ProteinMPNNResult from 'src/features/workflow/components/jobs/ProteinMpnnResult.vue';
import {
  nolabs__application__proteinmpnn__api_models__JobResponse,
  ProteinContentResponse,
  GetJobState,
  nolabs__application__proteinmpnn__api_models__SetupJobRequest
} from 'src/refinedApi/client';
import {
  getProteinMPNNJobApi,
  getProteinContent,
  getJobStatus,
  changeJobName,
  setupProteinMPNNJob,
  startProteinMPNNJob
} from 'src/features/workflow/refinedApi';

export default defineComponent({
  name: 'ProteinMPNNJob',
  props: {
    jobId: String,
  },
  data() {
    return {
      experimentId: null as string | null,
      job: null as nolabs__application__proteinmpnn__api_models__JobResponse | null,
      protein: null as ProteinContentResponse | null,
      jobStatus: null as GetJobState | null,
      editableJobName: '' as string,
      samplingTemp: null as number | null,
      numSeqPerTarget: null as number | null,
      isHomomer: false as boolean
    };
  },
  computed: {
    jobHasGeneratedData(): boolean | null {
      return this.job && this.job.result && this.job.result.length > 0;
    },
  },
  async mounted() {
    this.experimentId = this.$route.params.experimentId as string;

    this.job = await getProteinMPNNJobApi(this.jobId as string);
    if (this.job) {
      this.editableJobName = this.job.job_name || '';
      this.samplingTemp = this.job.sampling_temp || null;
      this.numSeqPerTarget = this.job.num_seq_per_target || null;
      this.isHomomer = this.job.is_homomer || false;
    }

    this.protein = await getProteinContent(this.job?.protein_id);

    this.jobStatus = await getJobStatus(this.jobId as string);
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
    async updateSamplingTemp() {
      if (this.job && this.samplingTemp !== this.job.sampling_temp) {
        const request: nolabs__application__proteinmpnn__api_models__SetupJobRequest = {
          experiment_id: this.experimentId as string,
          protein_id: this.job.protein_id,
          sampling_temp: this.samplingTemp,
          num_seq_per_target: this.numSeqPerTarget,
          is_homomer: this.isHomomer,
          job_id: this.job.job_id,
          job_name: this.job.job_name,
        };
        try {
          await setupProteinMPNNJob(request);
          this.$q.notify({
            type: 'positive',
            message: 'Sampling temperature updated successfully.',
          });
        } catch (error) {
          this.$q.notify({
            type: 'negative',
            message: 'Failed to update sampling temperature.',
          });
        }
      }
    },
    async updateNumSeqPerTarget() {
      if (this.job && this.numSeqPerTarget !== this.job.num_seq_per_target) {
        const request: nolabs__application__proteinmpnn__api_models__SetupJobRequest = {
          experiment_id: this.experimentId as string,
          protein_id: this.job.protein_id,
          sampling_temp: this.samplingTemp,
          num_seq_per_target: this.numSeqPerTarget,
          is_homomer: this.isHomomer,
          job_id: this.job.job_id,
          job_name: this.job.job_name,
        };
        try {
          await setupProteinMPNNJob(request);
          this.$q.notify({
            type: 'positive',
            message: 'Number of sequences per target updated successfully.',
          });
        } catch (error) {
          this.$q.notify({
            type: 'negative',
            message: 'Failed to update number of sequences per target.',
          });
        }
      }
    },
    async updateIsHomomer() {
      if (this.job && this.isHomomer !== this.job.is_homomer) {
        const request: nolabs__application__proteinmpnn__api_models__SetupJobRequest = {
          experiment_id: this.experimentId as string,
          protein_id: this.job.protein_id,
          sampling_temp: this.samplingTemp,
          num_seq_per_target: this.numSeqPerTarget,
          is_homomer: this.isHomomer,
          job_id: this.job.job_id,
          job_name: this.job.job_name,
        };
        try {
          await setupProteinMPNNJob(request);
          this.$q.notify({
            type: 'positive',
            message: 'Is Homomer updated successfully.',
          });
        } catch (error) {
          this.$q.notify({
            type: 'negative',
            message: 'Failed to update Is Homomer.',
          });
        }
      }
    },
    async startJob() {
      try {
        await startProteinMPNNJob(this.jobId as string);
        this.$q.notify({
          type: 'positive',
          message: 'Job started successfully.',
        });
        // Refresh the job status
        this.jobStatus = await getJobStatus(this.jobId as string);
      } catch (error) {
        this.$q.notify({
          type: 'negative',
          message: 'Failed to start the job.',
        });
      }
    }
  },
  components: {
    ProteinMPNNResult,
    QBtn,
    QToggle,
  },
});
</script>
