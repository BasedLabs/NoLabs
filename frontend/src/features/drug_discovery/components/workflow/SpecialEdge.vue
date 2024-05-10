<template>
  <svg class="special-edge" :style="edgeStyle">
    <path :d="edgePath" fill="none" stroke="#ffffff" stroke-width="2"></path> <!-- Change stroke color to white -->
  </svg>
</template>

<script>
export default {
  props: {
    edge: {
      type: Object,
      required: true
    }
  },
  computed: {
    edgePath() {
      // Generate path based on edge data if needed
      return `M${this.edge.sourcePosition.x},${this.edge.sourcePosition.y} L${this.edge.targetPosition.x},${this.edge.targetPosition.y}`;
    },
    edgeStyle() {
      // You can customize edge style here
      return {
        position: 'absolute',
        pointerEvents: 'none',
        zIndex: -1
      };
    }
  }
};
</script>

<style scoped>
.special-edge {
  width: 100%;
  height: 100%;
}
</style>
