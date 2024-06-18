<template>
  <q-card dark bordered>
    <q-card-section>
      <LigandsListNodeContent v-if="experimentId" :experiment-id="experimentId" :nodeId="nodeId" />
    </q-card-section>
    <q-card-section>
      <NodeHandles :nodeId="nodeId" :outputs="useNodesData?.data.outputs" />
    </q-card-section>
    <q-card-actions align="around">
      <q-btn @click="openDialogue" label="Extended view"/>
      <q-btn @click="deleteNode" color="negative" icon="delete" label="Delete" />
    </q-card-actions>
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from "vue";
import { useRoute } from "vue-router";
import { Position } from '@vue-flow/core'
import NodeHandles from '../../NodeHandles.vue'
import LigandsListNodeContent from "./LigandListNodeContent.vue";
import { useWorkflowStore } from "src/features/drug_discovery/components/workflow/storage";
import ProteinListNodeContent from "../proteins/ProteinListNodeContent.vue";

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
    nodeId: {
      type: String,
      required: true,
    },
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
    openDialogue() {
      this.onOpenDialog(this.nodeId);
    }
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
