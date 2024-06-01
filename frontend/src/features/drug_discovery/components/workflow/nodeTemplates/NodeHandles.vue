<template>
    <div class="q-mt-md">
      <q-item-label class="text-bold">Inputs:</q-item-label>
      <div v-for="(input, index) in inputs" :key="'input-' + index" class="row no-wrap items-center q-my-xs">
        <Handle type="target" :position="Position.Left" :id="`${nodeId}-input-${input}`" :style="{ position: 'relative', left: '50%', zIndex: 10 }" />
        <q-item-label class="q-ml-xs">{{ input }}</q-item-label>
      </div>
  
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
  .q-my-xs {
    margin-top: 4px;
    margin-bottom: 4px;
  }
  </style>
  