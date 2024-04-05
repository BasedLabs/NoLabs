<script lang="ts">
import {defineComponent, PropType} from 'vue'

export default defineComponent({
  name: "ExperimentControlButtons",
  props: {
    experiment: {
      type: Object as PropType<Experiment>,
      required: true
    },
    startExperiment: {
      type: Function as PropType<() => Promise<void>>,
      required: true
    },
    stopExperiment: {
      type: Function as PropType<() => Promise<void>>,
      required: true
    },
    resumeExperiment: {
      type: Function as PropType<() => Promise<void>>,
      required: true
    }
  },
  computed: {
    showRunButton() {
      return !this.experiment.running;
    },
    showStopButton() {
      return this.experiment.running;
    },
    showResumeButton() {
      return this.experiment.learningCompleted && !this.experiment.running;
    }
  }
})
</script>

<template>
  <q-btn align="between" outline size="lg" class="q-mx-md" color="positive" v-if="showRunButton" icon="biotech" @click="startExperiment">Start</q-btn>
  <q-btn align="between" outline size="lg" class="q-mx-md" color="negative" v-if="showStopButton" square dense @click="stopExperiment" icon="stop">Stop experiment</q-btn>
  <q-btn align="between" outline size="lg" class="q-mx-md" color="info" v-if="showResumeButton" square dense @click="resumeExperiment">
    Generate more smiles
  </q-btn>
</template>
