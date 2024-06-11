<template>
  <q-card class="q-pa-md" style="width: 400px">
    <div class="row no-wrap items-center">
      <q-btn @click="deleteNode" class="q-absolute-top-right" flat icon="delete" color="negative" />
      <q-space />
      <q-btn v-if="nodeData?.data.error" color="red" class="q-pm-md">
        <q-icon left name="warning" />
        <div>{{ nodeData?.data.error }}</div>
      </q-btn>
      <q-btn v-if="false" icon="settings" @click="openSettings"> </q-btn>
    </div>

    <JobNodeContent :nodeId="nodeId" :name="nodeData?.name" :description="nodeData?.description"/>

    <NodeHandles :nodeId="nodeId" :inputs="nodeData?.data.inputs" :outputs="nodeData?.data.outputs" />

    
    <q-btn @click="openDialogue" class="full-width q-pa-md" icon="open_in_new" label="Extended view"> </q-btn>
  </q-card>
</template>

<script>
import JobNodeContent from './JobNodeContent.vue';
import NodeHandles from './NodeHandles.vue';
import { useWorkflowStore } from 'src/features/drug_discovery/components/workflow/storage';

export default {
  name: 'JobNode',
  components: {
    JobNodeContent,
    NodeHandles
  },
  props: {
    nodeId: {
      type: String,
      required: true
    },
    onDeleteNode: Function,
    onOpenSettings: Function,
    onOpenDialog: Function
  },
  data() {
    return {
      nodeData: null
    };
  },
  created() {
    const workflowStore = useWorkflowStore();
    this.nodeData = workflowStore.getNodeById(this.nodeId);
  },
  methods: {
    deleteNode() {
      if (this.onDeleteNode) {
        this.onDeleteNode(this.nodeId);
      }
    },
    openSettings() {
      if (this.onOpenSettings) {
        this.onOpenSettings(this.nodeId);
      }
    },
    openDialogue() {
      if (this.onOpenDialog) {
        this.onOpenDialog(this.nodeId);
      }
    }
  }
};
</script>

<style scoped>
.job-section {
  border: 1px dashed #ccc;
  padding: 10px;
}
</style>
