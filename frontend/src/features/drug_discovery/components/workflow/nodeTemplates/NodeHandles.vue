<template>
    <div class="q-mt-md">
      <q-item-label class="text-bold q-mt-md q-mb-md">Inputs:</q-item-label>
      <div v-for="(input, index) in inputs" :key="'input-' + index" class="row no-wrap items-center q-pa-sm q-mb-sm q-border rounded-border q-shadow-2 bg-grey-3 input-output-tab">
        <Handle type="target" :position="Position.Left" :id="`${nodeId}-input-${input}`" class="large-handle" />
        <q-item-label class="q-mx-auto text-h6 text-black">{{ input }}</q-item-label>
      </div>
  
      <q-item-label class="text-bold q-mt-md">Outputs:</q-item-label>
      <div v-for="(output, index) in outputs" :key="'output-' + index" class="row no-wrap items-center q-pa-sm q-mt-md q-border rounded-border q-shadow-2 bg-grey-3 input-output-tab">
        <q-item-label class="q-mx-auto text-h6 text-black">{{ output }}</q-item-label>
        <Handle type="source" :position="Position.Right" :id="`${nodeId}-output-${output}`" class="large-handle" />
      </div>
    </div>
  </template>
  
  <script>
  import { defineComponent } from 'vue';
  import { Handle, Position } from "@vue-flow/core";
  
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
        default: () => []
      },
      outputs: {
        type: Array,
        default: () => []
      }
    },
    computed: {
      Position() {
        return Position;
      }
    }
  });
  </script>
  
  <style scoped>
  .text-h6 {
    font-size: 1.5rem; /* Adjust font size as needed */
  }
  
  .large-handle {
    width: 15px;
    height: 15px;
    position: relative;
    top: 5px;
    z-index: 10;
  }
  
  .input-output-tab {
    display: flex;
    height: 50px; /* Ensure a consistent height for centering */
  }
  </style>
  