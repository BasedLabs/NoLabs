<template>
  <q-card dark bordered>
    <q-card-section>
      <ProteinListNodeContent v-if="experimentId" :experiment-id="experimentId" :nodeId="nodeId" />
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
import { Notify, QSpinnerOrbit } from "quasar";
import { Position, useNodesData } from "@vue-flow/core";
import ProteinListNodeContent from "./ProteinListNodeContent.vue";
import NodeHandles from '../../NodeHandles.vue'
import { useWorkflowStore } from "src/features/drug_discovery/components/workflow/storage";
import JobNodeContent from "../../JobNodeContent.vue";

export default defineComponent({
  name: "ProteinListNode",
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
    ProteinListNodeContent,
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
    }
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
    this.$q.loading.show({
      spinner: QSpinnerOrbit,
      message: `Fetching targets for experiment`
    });
    try {
      //await store.fetchTargetsForExperiment(this.experimentId);
    } catch (e) {
      Notify.create({
        type: "negative",
        closeBtn: 'Close',
        message: 'Error fetching targets or ligands: ' + e
      });
    } finally {
      this.loading = false;
      this.$q.loading.hide();
    }
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
