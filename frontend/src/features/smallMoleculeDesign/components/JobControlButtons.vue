<script lang="ts">
import {defineComponent, PropType} from 'vue'

export default defineComponent({
  name: "JobControlButtons",
  props: {
    job: {
      type: Object as PropType<Job>,
      required: true
    },
    startJob: {
      type: Function as PropType<() => Promise<void>>,
      required: true
    },
    stopJob: {
      type: Function as PropType<() => Promise<void>>,
      required: true
    },
    startSampling: {
      type: Function as PropType<() => Promise<void>>,
      required: true
    },
    saveParameters: {
      type: Function as PropType<() => Promise<void>>,
      required: true
    }
  },
  computed: {
    showRunButton() {
      return !this.job.running;
    },
    showStopButton() {
      return this.job.running;
    },
    showSamplingButton() {
      return this.job.samplingAllowed && !this.job.running;
    }
  }
})
</script>

<template>
  <q-btn align="between" outline size="md" class="q-mx-md" color="positive" v-if="showRunButton" icon="biotech"
         @click="startJob">Start learning
    <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
      Start AI learning with provided parameters and protein target
    </q-tooltip>
  </q-btn>
  <q-btn align="between" outline size="md" class="q-mx-md" color="negative" v-if="showStopButton" square dense
         @click="stopJob" icon="stop">Stop job
  </q-btn>
  <q-btn align="between" outline size="md" class="q-mx-md" color="positive" v-if="showSamplingButton" square dense
         @click="startSampling">
    Start sampling
    <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
      Generate more molecules on trained AI
    </q-tooltip>
  </q-btn>
  <q-btn align="between" outline size="md" class="q-mx-md" color="positive" square dense
         @click="saveParameters">Save parameters
  </q-btn>
</template>
