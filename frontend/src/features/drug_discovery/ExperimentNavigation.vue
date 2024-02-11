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
        <q-item-label class="text-h5 q-pa-md">ExperimentId: {{ this.$route.params.experimentId }}</q-item-label>
        <div class="row no-wrap items-center q-mt-md text-white rounded-borders">
          <q-space />
          <q-btn flat label="Continue" @click="openNextStep('Upload ligands')" />
        </div>
        <router-view></router-view>
      </q-step>

      <q-step
        :name="2"
        title="Add ligands"
        :color="step > 2 ? 'info' : 'grey'"
        icon=hub
        :done="step > 2"
      >
        <q-item-label class="text-h5 q-pa-md">ExperimentId: {{ this.$route.params.experimentId }}</q-item-label>
        <div class="row no-wrap items-center q-mt-md text-white rounded-borders">
          <q-btn flat label="Back" @click="openPreviousStep('Upload targets')" />
          <q-space />
          <q-btn flat label="Continue" @click="openNextStep('Run docking')" />
        </div>
         <router-view></router-view>
      </q-step>

      <q-step
        :name="3"
        title="Run docking"
        :color="step > 3 ? 'info' : 'grey'"
        icon="gesture"
        :done="step > 3"
      >
        <q-item-label class="text-h5">ExperimentId: {{ this.$route.params.experimentId }}</q-item-label>
        <div class="row no-wrap items-center q-mt-md text-white rounded-borders">
          <q-btn flat label="Back" @click="openPreviousStep('Upload ligands')" />
        </div>
          <router-view></router-view>
      </q-step>
    </q-stepper>
  </q-page>
</template>

<script>
export default {
  name: "ExperimentNavigation",
  data() {
    return {
      step: 1, // This can be a reactive property based on the current route if needed
    };
  },
  mounted() {
    document.title = 'DTI';
    this.setStepBasedOnRoute();
  },
  methods: {
    setStepBasedOnRoute() {
      switch (this.$route.name) {
        case 'Upload targets':
          this.step = 1;
          break;
        case 'Upload ligands':
          this.step = 2;
          break;
        case 'Run docking':
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
    isStepDone(stepName) {
      // Implement logic to determine if a step is done based on your application's state
      // For example, checking if the current route's name matches the stepName
      return this.$route.name === stepName;
    }
  }
};
</script>
