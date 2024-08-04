<template>
  <div class="q-mt-md">
    <!-- Default job values section -->
    <div v-if="numberInputEntries.length" class="q-pa-sm q-mb-sm q-border rounded-border q-shadow-2 bg-grey-7">
      <div class="text-h6 text-white q-pl-sm q-mb-sm">Default job values</div>
      <div
        v-for="([key, input], index) in numberInputEntries"
        :key="'input-' + index"
        class="row no-wrap items-center q-pa-sm q-mb-sm input-output-tab"
      >
        <q-input
          v-model.number="currentInputValues[key]"
          @update:model-value="onInputChange(key, $event)"
          type="number"
          dense
          class="text-h6 text-white full-width-number-input"
          :label="input.title"
        />
      </div>
    </div>

    <!-- Node handles section -->
    <div
      v-for="([key, input], index) in nonNumberInputEntries"
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
import {defineComponent} from 'vue';
import {Handle, Position} from "@vue-flow/core";
import {useWorkflowStore} from '../storage';

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
  data() {
    return {
      currentInputValues: this.initializeInputValues()
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
    numberInputEntries() {
      return this.inputEntries.filter(([key, input]) => this.isNumberInput(input));
    },
    nonNumberInputEntries() {
      return this.inputEntries.filter(([key, input]) => !this.isNumberInput(input));
    }
  },
  methods: {
    isHandleConnected(handleId) {
      return this.edges.some(edge => edge.sourceHandle === handleId || edge.targetHandle === handleId);
    },
    isNumberInput(input) {
      return input.type === 'integer' || input.type === 'number';
    },
    initializeInputValues() {
      const values = {};
      Object.entries(this.inputs).forEach(([key, input]) => {
        if (this.isNumberInput(input)) {
          values[key] = input.default || 0;
        }
      });
      return values;
    },
    onInputChange(key, value) {
      const workflowStore = useWorkflowStore();
      this.currentInputValues[key] = value;
      workflowStore.setInputValue(this.nodeId, key, value);
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

.full-width-number-input .q-input-target {
  width: 100%;
  text-align: right;
}

.full-width-number-input .q-field {
  width: 100%;
}
</style>
