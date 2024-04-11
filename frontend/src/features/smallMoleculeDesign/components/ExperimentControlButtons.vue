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
    startSampling: {
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
    showSamplingButton() {
      return this.experiment.samplingAllowed && !this.experiment.running;
    }
  }
})
</script>

<template>
  <q-btn align="between" outline size="md" class="q-mx-md" color="positive" v-if="showRunButton" icon="biotech"
         @click="startExperiment">Start learning
    <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
      Start AI learning with provided parameters and protein target
    </q-tooltip>
  </q-btn>
  <q-btn align="between" outline size="md" class="q-mx-md" color="negative" v-if="showStopButton" square dense
         @click="stopExperiment" icon="stop">Stop experiment
  </q-btn>
  <q-btn align="between" outline size="md" class="q-mx-md" color="positive" v-if="showSamplingButton" square dense
         @click="startSampling">
    Start sampling
    <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
      Generate more molecules on trained AI
    </q-tooltip>
  </q-btn>
</template>
