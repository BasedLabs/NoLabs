<template>
  <q-card>
    <!-- Job Title with Icon -->
    <q-card-section>
      <div class="row items-center">
        <q-img
          src="/Adaptyvbio_small_logo.png"
          alt="Adaptyv Bio Logo"
          style="width: 50px; height: auto; margin-right: 10px;"
        />
        <div class="text-h6">
          AdaptyvBio Protein Affinity Characterization Job
        </div>
      </div>
    </q-card-section>

    <!-- Job Configuration -->
    <q-card-section>
      <q-list>
        <!-- Designs -->
        <q-item>
          <q-item-section>
            <q-item-label>Number of Designs</q-item-label>
          </q-item-section>
          <q-item-section>
            <div>
              <q-slider
                v-model="numDesigns"
                :min="24"
                :max="500"
                dense
                @change="fetchPriceEstimate"
                label-always
                markers
                :step="1"
                color="info"
              />
              <div class="row justify-between text-caption">
                <div>{{ 24 }}</div>
                <div>{{ 500 }}</div>
              </div>
            </div>
          </q-item-section>
        </q-item>
        <!-- Amino Acids -->
        <q-item>
          <q-item-section>
            <q-item-label>Amino Acids (AA)</q-item-label>
          </q-item-section>
          <q-item-section>
            <div>
              <q-slider
                v-model="numAA"
                :min="1"
                :max="300"
                dense
                @change="fetchPriceEstimate"
                label-always
                markers
                :step="1"
                color="info"
              />
              <div class="row justify-between text-caption">
                <div>{{ 1 }}</div>
                <div>{{ 300 }}</div>
              </div>
            </div>
          </q-item-section>
        </q-item>
        <!-- Replicates Per Design -->
        <q-item>
          <q-item-section>
            <q-item-label>Replicates per Design</q-item-label>
          </q-item-section>
          <q-item-section>
            <div>
              <q-slider
                v-model="replicatesPerDesign"
                :min="1"
                :max="5"
                dense
                @change="fetchPriceEstimate"
                label-always
                markers
                :step="1"
                color="info"
              />
              <div class="row justify-between text-caption">
                <div>{{ 1 }}</div>
                <div>{{ 5 }}</div>
              </div>
            </div>
          </q-item-section>
        </q-item>
        <!-- Email -->
        <q-item>
          <q-item-section>
            <q-item-label>Email for Quote</q-item-label>
          </q-item-section>
          <q-item-section>
            <q-input
              v-model="email"
              type="email"
              outlined
              dense
              clearable
              placeholder="example@gmail.com"
              @change="updateJob"
              :rules="[
                val => !!val || 'Email is required',
                val => /.+@.+\..+/.test(val) || 'Invalid email',
              ]"
            />
          </q-item-section>
        </q-item>
        <!-- Target -->
        <q-item>
          <q-item-section>
            <q-item-label>Target</q-item-label>
          </q-item-section>
          <q-item-section>
            <q-select
              v-model="selectedTarget"
              use-input
              :options="targetOptions"
              option-label="name"
              option-value="name"
              outlined
              dense
              @filter="onTargetFilter"
              :loading="isFetchingTargets"
              clearable
              input-debounce="1000"
              placeholder="Search target..."
              @update:model-value="fetchPriceEstimate"
            />
          </q-item-section>
        </q-item>
        <!-- Price Estimate -->
        <q-item>
          <q-item-section>
            <q-item-label>Price Estimate</q-item-label>
          </q-item-section>
          <q-item-section>
            <div v-if="isFetchingEstimate">
              <q-spinner size="24px" color="info" />
            </div>
            <div v-else>{{ priceEstimate }}</div>
          </q-item-section>
        </q-item>
      </q-list>
    </q-card-section>

    <!-- Submit Button -->
    <q-card-section v-if="!isJobSubmitted">
      <q-btn
        label="Submit for Review"
        color="primary"
        @click="submitJob"
        :disable="!readyForSubmission || isFetchingEstimate"
      />
    </q-card-section>

    <!-- Collaboration Phrase -->
    <q-card-section>
      <div class="text-caption text-center">
        For protein synthesis, we are collaborating with
        <a
          href="https://www.adaptyvbio.com/"
          target="_blank"
          class="text-info"
        >
          Adaptyv Bio
        </a>.
      </div>
    </q-card-section>
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import {
  QBtn,
  QInput,
  QSelect,
  QSlider,
  debounce,
  QImg,
  QSpinner,
} from 'quasar';
import {
  getAdaptyvBioEstimates,
  getProteinAffinityCharacterizationJob,
  listAvailableAdaptyvBioTargets,
  sendProteinAffinityCharacterizationRequest,
  setupProteinAffinityCharacterizationJob,
} from '../../refinedApi';

export type TargetResponse = {
  id: string;
  name: string;
  description: string;
  swissprot_id: string;
};

export type JobResponse = {
  job_id: string;
  job_name: string;
  number_of_designs: number;
  dna_length: number;
  replicates: number;
  report_email: string;
  target_id: string;
  swissprot_id: string;
  cart_total: number;
  session_url: string;
};

export default defineComponent({
  name: 'AdaptyvBioJob',
  data() {
    return {
      jobId: '', // Replace with the actual job ID
      numDesigns: 24,
      numAA: 178,
      replicatesPerDesign: 1,
      email: '',
      targetOptions: [] as TargetResponse[],
      selectedTarget: null as TargetResponse | null,
      priceEstimate: 'Loading...',
      estimateValue: 0, // Stores the numeric estimate
      isFetchingTargets: false,
      isFetchingEstimate: true, // Indicates if estimates are being fetched
      isJobSubmitted: false
    };
  },
  computed: {
    isEmailValid() {
      const email = this.email;
      if (!email) return false;
      const emailRegex = /.+@.+\..+/;
      return emailRegex.test(email);
    },
    readyForSubmission() {
      return (
        !this.isFetchingEstimate &&
        this.isEmailValid &&
        this.selectedTarget != null
      );
    },
  },
  methods: {
    async fetchJobParameters() {
      try {
        const job: JobResponse = await getProteinAffinityCharacterizationJob(
          this.jobId
        );

        // Set the component's state based on the fetched job parameters
        this.numDesigns = job.number_of_designs;
        this.numAA = job.dna_length;
        this.replicatesPerDesign = job.replicates;
        this.email = job.report_email;
        this.estimateValue = job.cart_total;
        this.priceEstimate = `$${job.cart_total} / estimate fetched`;
        this.selectedTarget = {
          id: job.target_id,
          name: job.target_id, // Replace with the target name if available
          description: '',
          swissprot_id: job.swissprot_id,
        };
        this.isJobSubmitted = job.submitted;
      } catch (error) {
        console.error('Failed to fetch job parameters:', error);
        this.$q.notify({
          type: 'negative',
          message: 'Failed to fetch job parameters.',
        });
      }
    },
    async fetchPriceEstimate() {
      this.isFetchingEstimate = true;
      try {
        const response = await getAdaptyvBioEstimates(this.jobId);
        this.priceEstimate = `$${response.total_price} / ${response.turnaround_time} weeks`;
        this.estimateValue = response.total_price;
        await this.updateJob();
      } catch (error) {
        this.priceEstimate = 'Failed to fetch estimate';
        console.error(error);
      } finally {
        this.isFetchingEstimate = false;
      }
    },
    async updateJob() {
      try {
        await setupProteinAffinityCharacterizationJob(
          this.jobId,
          this.numDesigns,
          this.numAA,
          this.replicatesPerDesign,
          this.email || null,
          this.selectedTarget?.id || null,
          this.estimateValue,
          this.selectedTarget?.swissprot_id || null
        );
      } catch (error) {
        this.$q.notify({
          type: 'negative',
          message: 'Failed to update the job.',
        });
        console.error(error);
      }
    },
    onTargetFilter(val, update, abort) {
      if (val === '') {
        update(() => {
          this.targetOptions = [];
        });
      } else {
        this.isFetchingTargets = true;
        this.debouncedFetchTargets(val, update);
      }
    },
    fetchTargets(val, update) {
      listAvailableAdaptyvBioTargets(val)
        .then((options: TargetResponse[]) => {
          update(() => {
            this.targetOptions = options;
          });
        })
        .catch((error) => {
          this.$q.notify({
            type: 'negative',
            message: 'Failed to fetch target options.',
          });
          console.error(error);
        })
        .finally(() => {
          this.isFetchingTargets = false;
        });
    },
    debouncedFetchTargets: debounce(function (val, update) {
      this.fetchTargets(val, update);
    }, 1000),
    async submitJob() {
      try {
        await sendProteinAffinityCharacterizationRequest(this.jobId);
        this.isJobSubmitted = true;
        this.$q.notify({
          type: 'positive',
          message: 'Job submitted successfully!',
        });
      } catch (error) {
        this.$q.notify({
          type: 'negative',
          message: 'Failed to submit the job.',
        });
        console.error(error);
      }
    },
  },
  async mounted() {
    // Fetch job parameters and then estimates
    this.jobId = this.$route.params.jobId as string;

    await this.fetchJobParameters(); // Fetch job details
    await this.fetchPriceEstimate(); // Fetch estimates based on the fetched job parameters
  },
  components: {
    QBtn,
    QSlider,
    QInput,
    QSelect,
    QImg,
    QSpinner,
  },
});
</script>
