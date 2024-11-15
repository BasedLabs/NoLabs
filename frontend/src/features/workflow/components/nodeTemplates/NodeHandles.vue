<template>
  <div class="q-mt-md">
    <!-- Default job values section -->
    <div
      v-if="defaultInputEntries.length"
      class="q-pa-sm q-mb-sm q-border rounded-border q-shadow-2 bg-grey-7"
    >
      <div class="text-h6 text-white q-pl-sm q-mb-sm">Default job values</div>
      <div
        v-for="([key, input], index) in defaultInputEntries"
        :key="'input-' + index"
        class="row no-wrap items-center q-pa-sm q-mb-sm input-output-tab"
      >
        <q-input
          v-model="currentInputValues[key]"
          @update:model-value="onInputChange(key, $event)"
          :type="input.type === 'string' ? 'text' : 'number'"
          hide-spin-buttons
          dense
          class="text-h6 text-white full-width-input"
          :label="input.title"
        />
      </div>
    </div>

    <!-- Node handles section -->
    <div
      v-for="([key, input], index) in otherInputEntries"
      :key="'input-' + index"
      class="row no-wrap items-center q-pa-sm q-mb-sm q-border rounded-border q-shadow-2 bg-grey-7 input-output-tab"
    >
      <Handle
        type="target"
        :position="Position.Left"
        :id="`${nodeId}-input-${key}`"
        :class="{
          'large-handle': true,
          'handle-connected': isHandleConnected(`${nodeId}-input-${key}`),
        }"
      />
      <q-item-label class="q-mx-auto text-h6 text-white">{{ input.title }}</q-item-label>
    </div>

    <!-- Output handles -->
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
        :class="{
          'large-handle': true,
          'handle-connected': isHandleConnected(`${nodeId}-output-${key}`),
        }"
      />
    </div>
  </div>
</template>

<script>
import { defineComponent } from 'vue';
import { Handle, Position } from '@vue-flow/core';
import { useWorkflowStore } from '../storage';

export default defineComponent({
  name: 'NodeHandles',
  components: {
    Handle,
  },
  props: {
    nodeId: {
      type: String,
      required: true,
    },
    inputs: {
      type: Object,
      default: () => ({}),
      required: false,
    },
    outputs: {
      type: Object,
      default: () => ({}),
      required: false,
    },
    defaults: {
      type: Array,
      default: () => [],
    },
  },
  data() {
    return {
      currentInputValues: this.initializeInputValues(),
    };
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
    },
    defaultInputEntries() {
      return this.inputEntries.filter(([key, input]) => this.isDefaultValueInput(input));
    },
    otherInputEntries() {
      return this.inputEntries.filter(([key, input]) => !this.isDefaultValueInput(input));
    },
  },
  methods: {
    isHandleConnected(handleId) {
      return this.edges.some(
        (edge) => edge.sourceHandle === handleId || edge.targetHandle === handleId
      );
    },
    isDefaultValueInput(input) {
      return ['integer', 'number', 'string'].includes(input.type);
    },
    initializeInputValues() {
      const values = {};
      Object.entries(this.inputs).forEach(([key, input]) => {
        if (this.isDefaultValueInput(input)) {
          values[key] =
            this.defaults.find((def) => def.target_path.includes(key))?.value ||
            input.default ||
            (input.type === 'string' ? '' : 0);
        }
      });
      return values;
    },
    onInputChange(key, value) {
      const workflowStore = useWorkflowStore();
      this.currentInputValues[key] = value;
      workflowStore.setInputValue(this.nodeId, key, value);
    },
  },
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
  align-items: center;
  height: 50px; /* Ensure a consistent height for centering */
  border-radius: 15px; /* Add border-radius for rounded corners */
}

.rounded-border {
  border-radius: 25px; /* Ensure this class has a border-radius as well */
}

.full-width-input {
  flex: 1;
}

.full-width-input .q-input-target {
  width: 100%;
}

.full-width-input .q-field {
  width: 100%;
}

/* Hide arrows in number input fields */
input[type='number']::-webkit-outer-spin-button,
input[type='number']::-webkit-inner-spin-button {
  -webkit-appearance: none;
  margin: 0;
}

input[type='number'] {
  -moz-appearance: textfield;
}
</style>
