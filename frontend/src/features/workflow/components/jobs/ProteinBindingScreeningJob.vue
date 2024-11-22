<template>
  <q-card>
    <q-card-section>
      <div class="text-h6">AdaptyvBio Job Configuration</div>
    </q-card-section>
    <q-card-section>
      <q-list>
        <!-- Experiment Type -->
        <q-item>
          <q-item-section>
            <q-item-label>Choose Experiment Type</q-item-label>
          </q-item-section>
          <q-item-section>
            <q-select
              v-model="experimentType"
              :options="experimentTypes"
              outlined
              dense
            />
          </q-item-section>
        </q-item>
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
      </q-list>
    </q-card-section>
    <q-card-section>
      <q-btn label="Submit for Review" color="primary" @click="submitJob" />
    </q-card-section>
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import { QBtn, QSelect, QSlider, QInput } from 'quasar';

export default defineComponent({
  name: 'AdaptyvBioJob',
  data() {
    return {
      experimentType: 'Protein Affinity Characterization',
      experimentTypes: [
        'Protein Affinity Characterization',
        'Enzyme Activity',
        'Structural Analysis',
        'Other',
      ],
      numDesigns: 1,
      numAA: 178,
      replicatesPerDesign: 1,
      email: '',
    };
  },
  methods: {
    submitJob() {
      if (!this.email) {
        this.$q.notify({
          type: 'negative',
          message: 'Email is required to submit the job.',
        });
        return;
      }

      // Here you can send the data to your API
      const jobData = {
        experimentType: this.experimentType,
        numDesigns: this.numDesigns,
        numAA: this.numAA,
        replicatesPerDesign: this.replicatesPerDesign,
        email: this.email,
      };

      console.log('Submitting job data:', jobData);

      // Simulate success
      this.$q.notify({
        type: 'positive',
        message: 'Job submitted successfully!',
      });
    },
  },
  components: {
    QBtn,
    QSelect,
    QSlider,
    QInput,
  },
});
</script>
