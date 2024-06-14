<template>
  <div class="q-mt-md">
    <q-item-label v-if="inputEntries.length > 0" class="text-bold q-mt-md q-mb-md text-h6">Inputs:</q-item-label>
    <div
      v-for="([key, input], index) in inputEntries"
      :key="'input-' + index"
      class="row no-wrap items-center q-pa-sm q-mb-sm q-border rounded-border q-shadow-2 bg-grey-7 input-output-tab"
    >
      <Handle
        type="target"
        :position="Position.Left"
        :id="`${nodeId}-input-${key}`"
        :class="{'large-handle': true, 'handle-connected': isHandleConnected(`${nodeId}-input-${key}`)}"
      />
      <q-item-label class="q-mx-auto text-h6 text-white">{{ input.title }}</q-item-label>
    </div>

    <q-item-label v-if="outputEntries.length > 0" class="text-bold q-mt-md text-h6">Outputs:</q-item-label>
    <div
      v-for="([key, output], index) in outputEntries"
      :key="'output-' + index"
      class="row no-wrap items-center q-pa-sm q-mt-md q-border bg-grey-7 rounded-border q-shadow-2 input-output-tab"
    >
      <q-item-label class="q-mx-auto text-h6 text-white">{{ output.title }}</q-item-label>
      <Handle
        type="source"
        :position="Position.Right"
        :id="`${nodeId}-output-${key}`"
        :class="{'large-handle': true, 'handle-connected': isHandleConnected(`${nodeId}-output-${key}`)}"
      />
    </div>
  </div>
</template>

<script>
import { defineComponent } from 'vue';
import { Handle, Position } from "@vue-flow/core";
import { useWorkflowStore } from 'src/features/drug_discovery/components/workflow/storage';

export default defineComponent({
  name: 'NodeHandles',
  components: {
    Handle
  },
  props: {
    nodeId: {
      type: String,
      required: true
    },
    inputs: {
      type: Object,
      default: () => ({}),
      required: false
    },
    outputs: {
      type: Object,
      default: () => ({}),
      required: false
    }
  },
  computed: {
    Position() {
      return Position;
    },
    edges() {
      const workflowStore = useWorkflowStore();
      return workflowStore.elements.edges;
    },
    inputEntries() {
      return Object.entries(this.inputs);
    },
    outputEntries() {
      return Object.entries(this.outputs);
    }
  },
  methods: {
    isHandleConnected(handleId) {
      return this.edges.some(edge => edge.sourceHandle === handleId || edge.targetHandle === handleId);
    }
  }
});
</script>

<style scoped>
.text-h6 {
  font-size: 1.5rem; 
}

.text-h5 {
  font-size: 1.75rem; 
}

.large-handle {
  width: 15px;
  height: 15px;
  position: relative;
  top: 6px;
  z-index: 10;
}

.handle-connected {
  background-color: rgb(0, 200, 255); /* Change this to your desired highlight color */
}

.input-output-tab {
  display: flex;
  height: 50px; /* Ensure a consistent height for centering */
  border-radius: 15px; /* Add border-radius for rounded corners */
}

.rounded-border {
  border-radius: 25px; /* Ensure this class has a border-radius as well */
}
</style>
