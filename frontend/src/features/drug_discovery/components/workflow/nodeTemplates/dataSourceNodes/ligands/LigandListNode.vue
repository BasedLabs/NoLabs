<template>
  <q-card class="q-pa-md" style="width: 400px">
    <div class="row no-wrap items-center">
      <q-btn @click="deleteNode" class="q-absolute-top-right" flat icon="delete" color="negative" />
      <q-space />
      <q-btn v-if="false" icon="settings" @click="openSettings"> </q-btn>
    </div>

    <LigandsListNodeContent v-if="experimentId" :experiment-id="experimentId" />
    <q-btn class="full-width" icon="open_in_new" label="Open in a new tab"> </q-btn>

    <NodeHandles :nodeId="nodeId" :outputs="useNodesData?.data.outputs" />
  </q-card>

</template>

<script lang="ts">
import { defineComponent } from "vue";
import { useRoute } from "vue-router";
import { Position } from '@vue-flow/core'
import NodeHandles from '../../NodeHandles.vue'
import LigandsListNodeContent from "./LigandListNodeContent.vue";
import { useWorkflowStore } from "src/features/drug_discovery/components/workflow/storage";

export default defineComponent({
  name: "LigandsListNode",
  computed: {
    Position() {
      return Position
    },
    useNodesData() {
      const workflowStore = useWorkflowStore();
      const nodeData = workflowStore.getNodeById(this.nodeId as string);
      return nodeData;
    }
  },
  components: {
    LigandsListNodeContent,
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
    onDeleteNode: {
      type: Function,
      required: true
    },
    onOpenSettings: {
      type: Function,
      required: true
    },
    onOpenDialog: {
      type: Function,
      required: true
    },
  },
  data() {
    return {
      loading: true,
      experimentId: null as string | null,
    };
  },
  async mounted() {
    const route = useRoute();

    this.experimentId = route.params.experimentId as string;
  },
  methods: {
    deleteNode() {
      this.onDeleteNode(this.nodeId);
    },
    openSettings() {
      this.onOpenSettings(this.nodeId);
    },
  }
})
</script>

<style scoped>
.special-node {
  border: 2px solid #ff6600;
  background-color: black;
  padding: 8px;
  border-radius: 4px;
  height: 400px;
  width: 400px
}

.special-node.input {
  border-color: #00bcd4;
  /* Change border color for input nodes */
}

.special-node.output {
  border-color: #4caf50;
  /* Change border color for output nodes */
}

.content {
  height: 100%;
  width: 100%;
}
</style>
