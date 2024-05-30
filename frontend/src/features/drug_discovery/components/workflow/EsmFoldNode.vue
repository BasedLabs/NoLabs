<template>
  <q-card class="q-pa-md" style="width: 400px">
    <div class="row no-wrap items-center">
      <q-btn @click="deleteNode" class="q-absolute-top-right" flat icon="delete" color="negative" />
      <q-space />
      <q-btn icon="settings" @click="openSettings"> </q-btn>
    </div>

    <EsmFoldNodeContent :results="results" />
    <q-btn @click="openDialogue" class="full-width" icon="open_in_new" label="Extended view"> </q-btn>
  </q-card>
  <div class="q-mt-md">
      <q-item-label class="text-bold">Inputs:</q-item-label>
      <div v-for="(input, index) in inputs" :key="'input-' + index" class="row no-wrap items-center q-my-xs">
        <Handle type="target" :position="Position.Left" :id="`${nodeId}-input-${input}`" :style="{ position: 'relative', left: '50%', zIndex: 10 }" />
        <q-item-label class="q-ml-xs">{{ input }}</q-item-label>
      </div>
    </div>

    <div class="q-mt-md">
      <q-item-label class="text-bold">Outputs:</q-item-label>
      <div v-for="(output, index) in outputs" :key="'output-' + index" class="row no-wrap items-center q-my-xs">
        <q-item-label class="q-mr-xs">{{ output }}</q-item-label>
        <Handle type="source" :position="Position.Right" :id="`${nodeId}-output-${output}`" :style="{ position: 'relative', left: '50%', zIndex: 10 }" />
      </div>
    </div>
</template>

<script>
import { defineComponent } from 'vue';
import { Handle, Position } from "@vue-flow/core";
import EsmFoldNodeContent from './EsmFoldNodeContent.vue';

export default defineComponent({
  name: 'EsmFoldNode',
  components: {
    EsmFoldNodeContent,
    Handle
  },
  props: {
    nodeId: String,
    inputs: {
      type: Array,
      default: () => []
    },
    outputs: {
      type: Array,
      default: () => []
    },
    onDeleteNode: Function,
    onOpenSettings: Function,
    onOpenDialog: Function,
    results: {
      type: Array,
      default: () => []
    },
  },
  computed: {
    Position() {
      return Position;
    }
  },
  methods: {
    deleteNode() {
      this.onDeleteNode(this.nodeId);
    },
    openSettings() {
      this.onOpenSettings(this.nodeId);
    },
    openDialogue() {
      this.onOpenDialog(this.nodeId);
    }
  }
});
</script>

<style scoped>
.job-section {
  border: 1px dashed #ccc;
  padding: 10px;
}
</style>
