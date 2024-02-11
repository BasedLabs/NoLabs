<template>
  <q-page class="q-pa-md ">
    <q-stepper
      v-model="step"
      ref="stepper"
      color="positive"
      animated
    >
      <q-step
        :name="1"
        title="Upload targets"
        icon=polymer
        color="info"
        :done="step > 1"
      >
        <ExperimentHeader v-if="experiment.metadata" :meta-data="experiment.metadata" :on-experiment-name-change-submit="onExperimentNameChange" />
        <div class="row no-wrap items-center q-mt-sm text-white rounded-borders">
          <q-space />
          <q-btn flat label="Continue" @click="openNextStep('UploadLigands')" />
        </div>
        <router-view></router-view>
      </q-step>

      <q-step
        :name="2"
        title="Add ligands"
        color="info"
        icon=hub
        :done="step > 2"
      >

        <div class="row no-wrap items-center q-mt-sm text-white rounded-borders">
          <q-btn flat label="Back" @click="openPreviousStep('UploadTargets')" />
          <q-space />
          <q-btn flat label="Continue" @click="openNextStep('RunDocking')" />
        </div>
         <router-view></router-view>
      </q-step>

      <q-step
        :name="3"
        title="Run docking"
        color="info"
        icon="gesture"
        :done="step > 3"
      >
        <div class="row no-wrap items-center q-mt-sm text-white rounded-borders">
          <q-btn flat label="Back" @click="openPreviousStep('UploadLigands')" />
        </div>
          <router-view></router-view>
      </q-step>
    </q-stepper>
  </q-page>
</template>

<script>
import ExperimentHeader from "src/features/drug_discovery/components/ExperimentHeader.vue";
import {useDrugDiscoveryStore} from "./storage";

export default {
  name: "ExperimentNavigation",
  components: {ExperimentHeader},
  data() {
    return {
      step: 1, // This can be a reactive property based on the current route if needed
      experiment: {},
    };
  },
  async mounted() {
    const store = useDrugDiscoveryStore();
    this.setStepBasedOnRoute();
    this.experiment.experimentId = this.$route.params.experimentId;
    this.experiment.metadata = await store.getExperimentMetaData(this.experiment.experimentId);
  },
  methods: {
    setStepBasedOnRoute() {
      switch (this.$route.name) {
        case 'UploadTargets':
          this.step = 1;
          break;
        case 'UploadLigands':
          this.step = 2;
          break;
        case 'RunDocking':
          this.step = 3;
          break;
        default:
          this.step = 1; // Default step if the route name doesn't match
      }
    },
    openNextStep(stepName) {
      this.$refs.stepper.next();
      this.$router.push({ name: stepName, params: { experimentId: this.$route.params.experimentId } });
    },
    openPreviousStep(stepName) {
      this.$refs.stepper.previous();
      this.$router.push({ name: stepName, params: { experimentId: this.$route.params.experimentId } });
    },
    async onExperimentNameChange(newExperimentName) {
      const store = useDrugDiscoveryStore();
      await store.changeExperimentName(this.experiment.experimentId, newExperimentName);
      this.experiment.name = newExperimentName;
    },
    isStepDone(stepName) {
      // Implement logic to determine if a step is done based on your application's state
      // For example, checking if the current route's name matches the stepName
      return this.$route.name === stepName;
    }
  }
};
</script>
