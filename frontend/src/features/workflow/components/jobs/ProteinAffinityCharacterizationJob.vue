<template>
  <q-card>
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
                :min="1"
                :max="4"
                dense
                @change="updateJob"
                label-always
                markers
                :step="1"
                color="info"
              />
              <div class="row justify-between text-caption">
                <div>{{ 1 }}</div>
                <div>{{ 4 }}</div>
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
                @change="updateJob"
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
                @change="updateJob"
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
              option-value="id"
              emit-value
              map-options
              outlined
              dense
              @filter="onTargetFilter"
              :loading="isFetchingTargets"
              clearable
              input-debounce="1000"
              placeholder="Search target..."
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
    <q-card-section>
      <div class="text-caption text-center">
        For protein synthesis, we are collaborating with
        <a href="https://www.adaptyvbio.com/" target="_blank">Adaptyv Bio</a>.
      </div>
    </q-card-section>
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import { QBtn, QInput, QSelect, QSlider, debounce } from 'quasar';
import {
  getAdaptyvBioEstimates,
  listAvailableAdaptyvBioTargets,
  sendProteinAffinityCharacterizationRequest,
  setupProteinAffinityCharacterizationJob,
} from '../../refinedApi';

export default defineComponent({
  name: 'AdaptyvBioJob',
  data() {
    return {
      jobId: 'your-job-id', // Replace with the actual job ID
      numDesigns: 1,
      numAA: 178,
      replicatesPerDesign: 1,
      email: '',
      targetOptions: [],
      selectedTarget: null,
      priceEstimate: 'Loading...',
      readyForSubmission: false,
      isFetchingTargets: false,
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
        await setupProteinAffinityCharacterizationJob(
          this.jobId,
          this.numDesigns,
          this.numAA,
          this.replicatesPerDesign,
          this.email,
          this.selectedTarget,
          0,
          '',
          ''
        );
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
        .then((options) => {
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
