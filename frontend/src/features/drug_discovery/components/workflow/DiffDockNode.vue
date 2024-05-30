<template>
  <q-card class="q-pa-md" style="width: 400px">
    <div class="row no-wrap items-center">
      <q-btn @click="deleteNode" class="q-absolute-top-right" flat icon="delete" color="negative" />
      <q-space />
      <q-btn icon="settings" @click="openSettings"> </q-btn>
    </div>

    <DiffDockNodeContent :jobs="jobs" :results="results" />
    <q-btn class="full-width" icon="open_in_new" label="Open in a new tab"> </q-btn>
    
    <div class="q-mt-md">
      <q-item-label class="text-bold">Inputs:</q-item-label>
      <div v-for="(input, index) in inputs" :key="'input-' + index" class="row no-wrap items-center q-my-xs">
        <q-item-label class="q-ml-xs">{{ input }}</q-item-label>
        <Handle type="target" :position="Position.Bottom" :id="`${nodeId}-input-${input}`" :style="{ position: 'relative', left: '50%' }" />
      </div>
    </div>

    <div class="q-mt-md">
      <q-item-label class="text-bold">Outputs:</q-item-label>
      <div v-for="(output, index) in outputs" :key="'output-' + index" class="row no-wrap items-center q-my-xs">
        <q-item-label class="q-mr-xs">{{ output }}</q-item-label>
        <Handle type="source" :position="Position.Bottom" :id="`${nodeId}-output-${output}`" :style="{ position: 'relative', left: '50%' }" />
      </div>
    </div>
  </q-card>
</template>

<script>
import {defineComponent} from 'vue';
import {Handle, Position} from "@vue-flow/core";
import DiffDockNodeContent from './DiffDockNodeContent.vue';

export default defineComponent({
  name: 'DiffDockNode',
  components: {
    DiffDockNodeContent,
    Handle
  },
  computed: {
    Position() {
      return Position;
    }
  },
  props: {
    nodeId: String,
    onDeleteNode: Function,
    onOpenSettings: Function,
    inputs: {
      type: Array,
      default: () => []
    },
    outputs: {
      type: Array,
      default: () => []
    },
    jobs: {
      type: Array,
      default: () => []
    },
    results: {
      type: Array,
      default: () => []
    },
  },
  methods: {
    deleteNode() {
      this.onDeleteNode(this.nodeId);
    },
    openSettings() {
      this.onOpenSettings(this.nodeId);
    },
  }
});
</script>

<style scoped>
.job-section {
  border: 1px dashed #ccc;
  padding: 10px;
}
</style>
