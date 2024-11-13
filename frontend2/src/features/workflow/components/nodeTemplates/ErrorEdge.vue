<script setup>
import { BaseEdge, EdgeLabelRenderer, getBezierPath } from '@vue-flow/core'
import { computed } from 'vue'
import {useWorkflowStore} from "../storage";

const props = defineProps({
  id: {
    type: String,
    required: true,
  },
  sourceX: {
    type: Number,
    required: true,
  },
  sourceY: {
    type: Number,
    required: true,
  },
  targetX: {
    type: Number,
    required: true,
  },
  targetY: {
    type: Number,
    required: true,
  },
  sourcePosition: {
    type: String,
    required: true,
  },
  targetPosition: {
    type: String,
    required: true,
  },
  data: {
    type: Object,
    required: false,
  },
  markerEnd: {
    type: String,
    required: false,
  },
  style: {
    type: Object,
    required: false,
  },
})

const path = computed(() => getBezierPath(props))

const workflowStore = useWorkflowStore()

const removeEdge = (id) => {
  workflowStore.onEdgeRemove(id)
}
</script>

<script>
export default {
  inheritAttrs: false,
}
</script>

<template>
  <BaseEdge
    :id="id"
    :style="{...style, strokeWidth: '4px'}"
    :path="path[0]"
    :marker-end="markerEnd"
    :label="data.text"
    :label-x="path[1]"
    :label-y="path[2]"
    :label-style="{ fill: 'white' }"
    :label-show-bg="true"
    :label-bg-style="{ fill: 'red' }"
    :label-bg-padding="[2, 4]"
    :label-bg-border-radius="2"
  />
  <EdgeLabelRenderer>
    <div
      :style="{
        pointerEvents: 'all',
        position: 'absolute',
        transform: `translate(-50%, -50%) translate(${path[1]}px,${path[2]}px)`,
      }"
      class="nodrag nopan"
    >
      <q-btn color="info" class="edgebutton" @click="removeEdge(id)">×</q-btn>
    </div>
  </EdgeLabelRenderer>
</template>

<style scoped>
.edgebutton {
  color: white;
  border: none;
  border-radius: 50%;
  width: 20px;
  height: 20px;
  cursor: pointer;
}
</style>
