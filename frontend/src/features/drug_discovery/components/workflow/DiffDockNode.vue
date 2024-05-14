<template>
  <q-card class="q-pa-md" style="width: 400px">
    <div class="row no-wrap items-center">
      <q-btn @click="deleteNode" class="q-absolute-top-right" flat icon="delete" color="negative" />
      <q-space />
      <q-btn icon="settings" @click="openSettings"> </q-btn>
    </div>

    <DiffDockNodeContent :jobs="jobs" :results="results" />
    <q-btn class="full-width" icon="open_in_new" label="Open in a new tab"> </q-btn>
    <Handle type="source" :position="Position.Right"/>
    <Handle type="target" :position="Position.Left"/>
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
