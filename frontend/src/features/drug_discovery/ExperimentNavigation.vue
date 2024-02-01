<template>
  <q-page class="q-pa-md">
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
        <div class="row no-wrap items-center q-mt-md text-white rounded-borders">
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
        <div class="row no-wrap items-center q-mt-md text-white rounded-borders">
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
        <div class="row no-wrap items-center q-mt-md text-white rounded-borders">
          <q-btn flat label="Back" @click="openPreviousStep('UploadLigands')" />
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
  methods: {
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
