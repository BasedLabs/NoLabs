<template>
  <q-card class="q-pa-md" style="width: 400px">
    <div class="row no-wrap items-center">
      <q-btn @click="deleteNode" class="q-absolute-top-right" flat icon="delete" color="negative" />
      <q-space />
      <q-btn v-if="nodeData?.data.error" color="red" class="q-pm-md">
        <q-icon left name="warning" />
        <div>{{ nodeData?.data.error }}</div>
      </q-btn>
      <q-btn icon="settings" @click="openSettings"> </q-btn>
    </div>

    <JobNodeContent :name="nodeData?.name" :description="nodeData?.description" :jobIds="nodeData?.jobIds" :results="results" />
    <q-btn @click="openDialogue" class="full-width" icon="open_in_new" label="Extended view"> </q-btn>

    <NodeHandles :nodeId="nodeId" :inputs="nodeData?.data.inputs" :outputs="nodeData?.data.outputs" />
  </q-card>
</template>

<script>
import { defineComponent, computed } from 'vue';
import JobNodeContent from './JobNodeContent.vue';
import NodeHandles from './nodeTemplates/NodeHandles.vue';
import { useWorkflowStore } from 'src/features/drug_discovery/components/workflow/storage';

export default defineComponent({
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
  setup(props) {
    const workflowStore = useWorkflowStore();
    const nodeData = workflowStore.getNodeById(props.nodeId);
    
    const deleteNode = () => {
      if (props.onDeleteNode) {
        props.onDeleteNode(props.nodeId);
      }
    };

    const openSettings = () => {
      if (props.onOpenSettings) {
        props.onOpenSettings(props.nodeId);
      }
    };

    const openDialogue = () => {
      if (props.onOpenDialog) {
        props.onOpenDialog(props.nodeId);
      }
    };

    return {
      nodeData,
      deleteNode,
      openSettings,
      openDialogue
    };
  }
});
</script>

<style scoped>
.job-section {
  border: 1px dashed #ccc;
  padding: 10px;
}
</style>
