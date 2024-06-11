<template>
  <div class="q-mt-md">
    <q-item-label v-if="inputs.length > 0" class="text-bold q-mt-md q-mb-md text-h6">Inputs:</q-item-label>
    <div
      v-for="(input, index) in inputs"
      :key="'input-' + index"
      class="row no-wrap items-center q-pa-sm q-mb-sm q-border rounded-border q-shadow-2 bg-grey-7 input-output-tab"
    >
      <Handle
        type="target"
        :position="Position.Left"
        :id="`${nodeId}-input-${input}`"
        :class="{'large-handle': true, 'handle-connected': isHandleConnected(`${nodeId}-input-${input}`)}"
      />
      <q-item-label class="q-mx-auto text-h6 text-white">{{ input }}</q-item-label>
    </div>

    <q-item-label v-if="outputs.length > 0" class="text-bold q-mt-md text-h6">Outputs:</q-item-label>
    <div
      v-for="(output, index) in outputs"
      :key="'output-' + index"
      class="row no-wrap items-center q-pa-sm q-mt-md q-border bg-grey-7 rounded-border q-shadow-2 input-output-tab"
    >
      <q-item-label class="q-mx-auto text-h6 text-white">{{ output }}</q-item-label>
      <Handle
        type="source"
        :position="Position.Right"
        :id="`${nodeId}-output-${output}`"
        :class="{'large-handle': true, 'handle-connected': isHandleConnected(`${nodeId}-output-${output}`)}"
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
      type: Array,
      default: () => [],
      required: false
    },
    outputs: {
      type: Array,
      default: () => [],
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
  font-size: 1.5rem; /* Adjust font size as needed */
}

.text-h5 {
  font-size: 1.75rem; /* Adjust font size as needed */
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
