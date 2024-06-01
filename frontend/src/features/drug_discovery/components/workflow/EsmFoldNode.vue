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
  
  <NodeHandles :nodeId="nodeId" :inputs="inputs" :outputs="outputs" />
</template>

<script>
import { defineComponent } from 'vue';
import { Position } from "@vue-flow/core";
import EsmFoldNodeContent from './EsmFoldNodeContent.vue';
import NodeHandles from './nodeTemplates/NodeHandles.vue'

export default defineComponent({
  name: 'EsmFoldNode',
  components: {
    EsmFoldNodeContent,
    NodeHandles
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
