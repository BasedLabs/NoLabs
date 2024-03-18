<template>
  <q-page class="bg-black q-pl-md q-pr-md">
        <BioBuddyChat v-if="bioBuddyEnabled" :experiment-id=" this.$route.params.experimentId" />
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
            :color="step > 0 ? 'info' : 'grey'"
            :done="step > 1"
          >
            <ExperimentHeader v-if="experiment.metadata" :meta-data="experiment.metadata" :on-experiment-name-change-submit="onExperimentNameChange" />
            <div class="row no-wrap items-center q-mt-sm text-white rounded-borders">
              <q-space />
              <q-btn flat label="Continue" @click="openNextStep('Upload ligands')" />
            </div>
            <router-view></router-view>
          </q-step>

          <q-step
            :name="2"
            title="Add ligands"
            :color="step > 1 ? 'info' : 'grey'"
            icon=hub
            :done="step > 2"
          >

        <div class="row no-wrap items-center q-mt-sm text-white rounded-borders">
          <q-btn flat color="info" label="Back" @click="openPreviousStep('Upload targets')" />
          <q-space />
          <q-btn flat color="info" label="Continue" @click="openNextStep('Run docking')" />
        </div>
         <router-view></router-view>
      </q-step>

          <q-step
            :name="3"
            title="Run docking"
            :color="step > 2 ? 'info' : 'grey'"
            icon="gesture"
            :done="step > 3"
          >
            <div class="row no-wrap items-center q-mt-sm text-white rounded-borders">
              <q-btn flat label="Back" @click="openPreviousStep('Upload ligands')" />
            </div>
              <router-view></router-view>
          </q-step>
        </q-stepper>
  </q-page>
</template>

<script lang="ts">
import ExperimentHeader from "src/features/drug_discovery/components/ExperimentHeader.vue";
import BioBuddyChat from "src/features/biobuddy/BioBuddyChat.vue";
import {useDrugDiscoveryStore} from "./storage";
import {defineComponent} from "vue";
import {checkBioBuddyEnabled} from "../biobuddy/api";


export default defineComponent({
  name: "ExperimentNavigation",
  components: {BioBuddyChat, ExperimentHeader},
  data() {
    return {
      step: 1, // This can be a reactive property based on the current route if needed
      experiment: {
        experimentId: null as string | null,
        metadata: null
      },
      splitterModel: 20,
      bioBuddyEnabled: false,
    };
  },
  async mounted() {
    const store = useDrugDiscoveryStore();
    this.setStepBasedOnRoute();
    this.experiment.experimentId = this.$route.params.experimentId as string;
    this.experiment.metadata = await store.getExperimentMetaData(this.experiment.experimentId);
    try {
      const response = await checkBioBuddyEnabled();
      this.bioBuddyEnabled = response.enabled;
    } catch (error) {
      console.error('Error checking BioBuddy enabled status:', error);
      this.bioBuddyEnabled = false;
    }
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
    openNextStep(stepName: string) {
      this.$refs.stepper.next();
      this.$router.push({ name: stepName, params: { experimentId: this.$route.params.experimentId } });
    },
    openPreviousStep(stepName: string) {
      this.$refs.stepper.previous();
      this.$router.push({ name: stepName, params: { experimentId: this.$route.params.experimentId } });
    },
    async onExperimentNameChange(newExperimentName: string) {
      const store = useDrugDiscoveryStore();
      await store.changeExperimentName(this.experiment.experimentId as string, newExperimentName);
    }
  }
}
)
</script>

<style>

body {
  overflow-x: hidden;
}

</style>
