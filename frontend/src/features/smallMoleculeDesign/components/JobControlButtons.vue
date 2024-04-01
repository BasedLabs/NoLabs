<script lang="ts">
import {defineComponent, PropType} from 'vue'

export default defineComponent({
  name: "jobControlButtons",
  props: {
    job: {
      type: Object as PropType<Job>,
      required: false
    },
    startJob: {
      type: Function as PropType<(jobId: string) => Promise<void>>,
      required: true
    },
    stopJob: {
      type: Function as PropType<(jobId: string) => Promise<void>>,
      required: true
    },
    resumeJob: {
      type: Function as PropType<(jobId: string) => Promise<void>>,
      required: true
    }
  },
  computed: {
    showRunButton() {
      return !this.job?.running;
    },
    showStopButton(){
      return this.job?.running;
    },
    showResumeButton() {
      return this.job?.learningCompleted && !this.job.running;
    }
  }
})
</script>

<template>
  <q-td auto-width v-if="showRunButton">
    <q-btn size="md" color="info" square dense @click="startJob(job!.id)" icon="play_arrow">Start</q-btn>
  </q-td>
  <q-td auto-width v-if="showStopButton">
    <q-btn size="md" color="negative" square dense @click="stopJob(job!.id)" icon="stop">Stop</q-btn>
  </q-td>
  <q-td auto-width v-if="showResumeButton">
    <q-btn size="md" color="info" square dense @click="resumeJob(job!.id)" icon="play_arrow">Generate more smiles</q-btn>
  </q-td>
</template>
