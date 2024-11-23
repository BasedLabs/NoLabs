<template>
  <q-card>
    <q-card-section>
      <div class="text-h6">AdaptyvBio Job Configuration</div>
    </q-card-section>
    <q-card-section>
      <q-list>
        <!-- Designs -->
        <q-item>
          <q-item-section>
            <q-item-label>Number of Designs</q-item-label>
          </q-item-section>
          <q-item-section>
            <q-slider
              v-model="numDesigns"
              :min="1"
              :max="4"
              dense
              @change="updateJob"
            />
          </q-item-section>
        </q-item>
        <!-- Amino Acids -->
        <q-item>
          <q-item-section>
            <q-item-label>Amino Acids (AA)</q-item-label>
          </q-item-section>
          <q-item-section>
            <q-slider
              v-model="numAA"
              :min="1"
              :max="300"
              dense
              @change="updateJob"
            />
          </q-item-section>
        </q-item>
        <!-- Replicates Per Design -->
        <q-item>
          <q-item-section>
            <q-item-label>Replicates per Design</q-item-label>
          </q-item-section>
          <q-item-section>
            <q-slider
              v-model="replicatesPerDesign"
              :min="1"
              :max="5"
              dense
              @change="updateJob"
            />
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
            />
          </q-item-section>
        </q-item>
        <!-- Target -->
        <q-item>
          <q-item-section>
            <q-item-label>Target</q-item-label>
          </q-item-section>
          <q-item-section>
            <q-input
              v-model="targetSearch"
              @change="searchTargets"
              dense
              outlined
              placeholder="Search target..."
            />
            <q-select
              v-model="selectedTarget"
              :options="targetOptions"
              option-label="name"
              option-value="id"
              emit-value
              map-options
              outlined
              dense
            />
          </q-item-section>
        </q-item>
        <!-- Price Estimate -->
        <q-item>
          <q-item-section>
            <q-item-label>Price Estimate</q-item-label>
          </q-item-section>
          <q-item-section>
            <div>{{ priceEstimate }}</div>
          </q-item-section>
        </q-item>
      </q-list>
    </q-card-section>
    <q-card-section>
      <q-btn
        label="Submit for Review"
        color="primary"
        @click="submitJob"
        :disable="!readyForSubmission"
      />
    </q-card-section>
  </q-card>
</template>

<script lang="ts">
import {defineComponent} from 'vue';
import {QBtn, QInput, QSelect, QSlider} from 'quasar';
import {ProteinAffinityCharacterizationService} from 'src/refinedApi/client';
import {
  getAdaptyvBioEstimates,
  listAvailableAdaptyvBioTargets,
  sendProteinAffinityCharacterizationRequest,
  setupProteinAffinityCharacterizationJob
} from "../../refinedApi";

export default defineComponent({
  name: 'AdaptyvBioJob',
  data() {
    return {
      jobId: 'your-job-id', // Replace with the actual job ID
      numDesigns: 1,
      numAA: 178,
      replicatesPerDesign: 1,
      email: '',
      targetSearch: '',
      targetOptions: [],
      selectedTarget: null,
      priceEstimate: 'Loading...',
      readyForSubmission: false,
    };
  },
  methods: {
    async fetchPriceEstimate() {
      try {
        const response = await getAdaptyvBioEstimates(this.jobId);
        this.priceEstimate = `$${response.total_price} / ${response.turnaround_time} weeks`;
      } catch (error) {
        this.priceEstimate = 'Failed to fetch estimate';
        console.error(error);
      }
    },
    async updateJob() {
      if (!this.email || !this.selectedTarget) {
        this.$q.notify({
          type: 'negative',
          message: 'Email and target must be set to update the job.',
        });
        return;
      }
      try {
        await setupProteinAffinityCharacterizationJob(this.jobId,
          this.numDesigns,
          this.numAA,
          this.replicatesPerDesign,
          this.email,
          this.selectedTarget,
          0,
          '',
          ''
        )
        this.readyForSubmission = true;
        this.$q.notify({
          type: 'positive',
          message: 'Job updated successfully.',
        });
        await this.fetchPriceEstimate();
      } catch (error) {
        this.$q.notify({
          type: 'negative',
          message: 'Failed to update the job.',
        });
        console.error(error);
      }
    },
    async searchTargets() {
      try {
        this.targetOptions = await listAvailableAdaptyvBioTargets(this.targetSearch);
      } catch (error) {
        this.$q.notify({
          type: 'negative',
          message: 'Failed to fetch target options.',
        });
        console.error(error);
      }
    },
    async submitJob() {
      try {
        await sendProteinAffinityCharacterizationRequest(this.jobId);
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
    // Fetch the initial price estimate when the component mounts
    await this.fetchPriceEstimate();
  },
  components: {
    QBtn,
    QSlider,
    QInput,
    QSelect,
  },
});
</script>
